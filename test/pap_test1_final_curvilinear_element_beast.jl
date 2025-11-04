using StaticArrays
using LinearAlgebra, Statistics
using Test

using CompScienceMeshes
import CompScienceMeshes: chart, measure, jacobian, refnodes, neighborhood,
                          coordtype, dimension, universedimension, vertextype,
                          quadpoints

using BEAST
using FastGaussQuadrature
using SpecialFunctions
using Gmsh

using LinearAlgebra: norm

"""
This test sweeps `p = 1:5` and validates the curvilinear high-order integration pipeline.

It assembles the 2D Helmholtz single-layer operator on a circular boundary, 
verifying the following for each polynomial order:

1. **Assembly symmetry** — ensures complex-symmetric matrices from BEAST’s integral operators.
2. **Sauter–Schwab vs. product-rule consistency** — compares DoubleNumSauterQStrat and DoubleNumQStrat.
3. **Quadrature refinement stability** — checks convergence with increasing quadrature order.
4. **Rayleigh-quotient eigenvalue accuracy** — validates mode-by-mode consistency against analytic solutions.
5. **Geometric contact classification** — confirms Sauter–Schwab edge/vertex/product rule detection.
6. **Optional basis equivalence** — (if enabled) compares BEAST’s native D0 basis with the 
   Gmsh-node-based `lagrangec0_gmsh_exact()` implementation.

This serves as a full integration test for curvilinear element geometry, 
basis construction, and high-order quadrature consistency in BEAST.jl.
"""



quad_measure(ch, n::Int=5) = begin
    q = CompScienceMeshes.quadpoints(ch, n)
    w = q isa Tuple{AbstractVector,AbstractVector} ? q[2] : last.(q)
    sum(w)
end

function lagrangec0_gmsh_exact(mesh)
    T = coordtype(mesh)
    P = vertextype(mesh)
    S = BEAST.Shape{T}

    dofs = Dict{Int, Vector{S}}()

    first_face = first(mesh.faces)
    NFcell = length(first_face)

    for (c, cell) in pairs(mesh.faces)
        @assert length(cell) == NFcell "Mixed-order faces not supported"
        @inbounds for r in 1:NFcell
            g = cell[r]
            push!(get!(dofs, g, Vector{S}()), S(c, r, T(1)))
        end
    end

    gids = sort!(collect(keys(dofs)))
    fns  = [dofs[g] for g in gids]
    pos  = [mesh.vertices[g] for g in gids]

    porder = NFcell - 1
    return BEAST.LagrangeBasis{porder, 0, NFcell}(mesh, fns, pos)
end

"""
    circle_curvilinear(radius, porder; h = 2π*radius/64)

Return a 1D-in-2D CurvilinearMesh of a circle with polynomial order `porder`
(using Gmsh high-order line elements). `h` controls target edge length.
"""
function circle_curvilinear(radius::Real, porder::Integer; h::Real = 2π*radius/64)
    @assert porder ≥ 1 "porder must be ≥ 1"
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity", 0)

        gmsh.model.add("circle_p$(porder)")
        s = gmsh.model.occ.addDisk(0.0, 0.0, 0.0, radius, radius)
        gmsh.model.occ.synchronize()

        # sizing + high-order
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
        gmsh.option.setNumber("Mesh.ElementOrder", porder)
        gmsh.option.setNumber("Mesh.HighOrderOptimize", porder >= 2 ? 2 : 0)
        gmsh.option.setNumber("Mesh.SecondOrderLinear", porder >= 2 ? 0 : 1)

        # physicals
        gmsh.model.addPhysicalGroup(2, [s], 1)
        gmsh.model.setPhysicalName(2, 1, "Domain")
        curves = [t[2] for t in gmsh.model.getBoundary([(2, s)], false, false, false) if t[1] == 1]
        gmsh.model.addPhysicalGroup(1, curves, 2)
        gmsh.model.setPhysicalName(1, 2, "Boundary")

        gmsh.model.mesh.generate(2)

        # nodes on boundary physical
        nb = gmsh.model.mesh.getNodesForPhysicalGroup(1, 2)
        nodeTags   = nb[1]
        nodeCoords = nb[2]
        @assert length(nodeCoords) == 3*length(nodeTags) "unexpected coords size"

        tag2idx = Dict{Int,Int}(Int(nodeTags[i]) => i for i in eachindex(nodeTags))

        verts = Vector{SVector{2,Float64}}(undef, length(nodeTags))
        @inbounds for i in eachindex(nodeTags)
            x = nodeCoords[3(i-1)+1]; y = nodeCoords[3(i-1)+2]
            verts[i] = SVector{2,Float64}(x, y)
        end

        NF = porder + 1
        faces = SVector{NF,Int}[]
        for c in curves
            types, elemTags, elemNodeTags = gmsh.model.mesh.getElements(1, c)
            @inbounds for k in eachindex(types)
                tags  = elemTags[k]
                nodes = elemNodeTags[k]
                ne = length(tags)
                ne == 0 && continue

                nNode_block = div(length(nodes), ne)
                nNode_block == NF || continue

                @assert length(nodes) == ne * nNode_block
                for e in 1:ne
                    off = (e-1) * nNode_block
                    ids = ntuple(j -> tag2idx[Int(nodes[off + j])], NF)
                    push!(faces, SVector{NF,Int}(ids))
                end
            end
        end
        @assert !isempty(faces) "No boundary line elements with p=$(porder) found."

        return CurvilinearMesh(verts, faces, porder)
    finally
        gmsh.finalize()
    end
end

# ----------------------------------------------------------------------
# Shared helpers
# ----------------------------------------------------------------------
edgehits_local(chs, cht; tol::Real=1e-12) = begin
    @assert CompScienceMeshes.dimension(chs) == 1
    @assert CompScienceMeshes.dimension(cht) == 1
    s0, s1 = CompScienceMeshes.nodes(chs)
    t0, t1 = CompScienceMeshes.nodes(cht)
    d00 = norm(s0 - t0); d01 = norm(s0 - t1)
    d10 = norm(s1 - t0); d11 = norm(s1 - t1)
    if (d00 ≤ tol && d11 ≤ tol) || (d01 ≤ tol && d10 ≤ tol)
        2
    elseif d00 ≤ tol || d01 ≤ tol || d10 ≤ tol || d11 ≤ tol
        1
    else
        0
    end
end

const SSQ1D = BEAST.SauterSchwabQuadrature1D

# ----------------------------------------------------------------------
# Main sweep over p = 1:5
# ----------------------------------------------------------------------
R = 1.0
N = 64
hs = 2π*R/N

for p in 1:5
    @testset "Curvilinear circle — order p=$p" begin
        # mesh + basis
        m = circle_curvilinear(R, p; h=hs)
        X = lagrangecxd0(m)                # BEAST's own D0 basis from the curved mesh
        # alt exact-DOF basis (works for any order as implemented)
        #X_exact = lagrangec0_gmsh_exact(m)

        # ---------------- Assembly & symmetry ----------------
        S = Helmholtz2D.singlelayer(wavenumber=1.0)

        Zq  = assemble(S, X, X;
                       threading = BEAST.Threading{:single},
                       quadstrat = BEAST.DoubleNumQStrat(15,14))
        @test norm(Zq - transpose(Zq)) / max(norm(Zq), eps()) < 1e-12

        Zss = assemble(S, X, X;
                       threading = BEAST.Threading{:single},
                       quadstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30))
        @test norm(Zss - transpose(Zss)) / max(norm(Zss), eps()) < 1e-12

        δ = norm(0.5*(Zq + transpose(Zq)) - 0.5*(Zss + transpose(Zss))) / max(norm(Zss), eps())
        @test δ < 5e-3

        # ---------------- Order refinement (stability) ----------------
        sym(A) = 0.5 .* (A .+ transpose(A))
        qref = 20
        Zref = sym(assemble(S, X, X;
                            threading = BEAST.Threading{:single},
                            quadstrat = BEAST.DoubleNumSauterQstrat(qref,qref, qref,qref, qref,qref)))

        for q in (10, 15, 20)
            Z = sym(assemble(S, X, X;
                             threading = BEAST.Threading{:single},
                             quadstrat = BEAST.DoubleNumSauterQstrat(q,q, q,q, q,q)))
            δq = norm(Z - Zref) / max(norm(Zref), eps())
            @info "p=$p  Sauter–Schwab q=$q vs ref(q=$qref)" δq
            @test δq < 5e-3
        end

        # ---------------- Rayleigh quotient checks ----------------
        κ = 3.0
        Sκ = Helmholtz2D.singlelayer(wavenumber = κ)
        Zκ = sym(assemble(Sκ, X, X;
                          threading = BEAST.Threading{:single},
                          quadstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)))

        MD0 = Diagonal([measure(chart(m,i)) for i in cells(m)])

        λ_analytic(n) = -0.5im * π * R * besselj(n, κ*R) * hankelh2(n, κ*R)

        function d0_mode(n::Int; q=12)
            c = ComplexF64[]
            for i in cells(m)
                ch = chart(m, i)
                ξ, w = gausslegendre(q)
                s = 0.0 + 0.0im
                len = measure(ch)
                @inbounds for k in eachindex(ξ)
                    ζ  = (ξ[k] + 1)/2
                    x  = map(ch, ζ)
                    θ  = atan(x[2], x[1]) % (2π)
                    s += cis(n*θ) * jacobian(ch, ζ) * 0.5 * w[k]
                end
                push!(c, s/len)
            end
            c
        end

        for n in 0:6
            c = d0_mode(n; q=12)
            λ_RQ = (c' * Zκ * c) / (c' * (MD0 * c))
            λ_ex = λ_analytic(n)
            @info "p=$p  n=$n  |λ_RQ-λ_ex|/|λ_ex|" rel = abs(λ_RQ-λ_ex)/max(abs(λ_ex), eps())
            @test isapprox(λ_RQ, λ_ex; rtol=3e-2)
        end

        # ---------------- SS contact classification ----------------
        is = first(cells(m))
        it = is == last(cells(m)) ? first(cells(m)) : is + 1
        τ  = chart(m, is)
        σ  = chart(m, it)
        k  = is + 5; k > last(cells(m)) && (k -= length(cells(m)))
        ρ  = chart(m, k)

        @test edgehits_local(τ, τ) == 2
        @test edgehits_local(τ, σ) == 1
        @test edgehits_local(τ, ρ) == 0

        # Rule selection
        g  = BEAST.refspace(X)
        qs = BEAST.DoubleNumSauterQstrat(10,10, 10,10, 10,10)
        qd = BEAST.quaddata(S, g, g, [τ], [σ], qs)

        @test BEAST.quadrule(S, g, g, 1, τ, 1, τ, qd, qs) isa SSQ1D.CommonEdge
        @test BEAST.quadrule(S, g, g, 1, τ, 1, σ, qd, qs) isa SSQ1D.CommonVertex

        k = mod1(is + 5, length(cells(m)))
        if k == is || k == it || k == mod1(is-1, length(cells(m)))
            k = mod1(is + length(cells(m)) ÷ 2, length(cells(m)))
        end
        ρ = chart(m, k)
        r = BEAST.quadrule(S, g, g, 1, τ, 1, ρ, qd, qs)
        @test !(r isa SSQ1D.CommonEdge) && !(r isa SSQ1D.CommonVertex)
    end
end
