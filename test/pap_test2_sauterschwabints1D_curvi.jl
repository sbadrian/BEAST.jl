using Test
using LinearAlgebra
using StaticArrays
using Gmsh
using BEAST
using CompScienceMeshes
import CompScienceMeshes: chart, measure, jacobian, coordtype, vertextype

"""
It sweeps p=1:5, assembles the 2D Helmholtz single-layer operator, checks symmetry, 
Sauter–Schwab vs product-rule consistency, optional basis equivalence, 
contact classification, and quadrature refinement.
"""


# ---------------------------
# Helpers you wanted to use
# ---------------------------
function lagrangec0_gmsh_exact(mesh)
    T = coordtype(mesh)
    P = vertextype(mesh)
    S = BEAST.Shape{T}

    dofs = Dict{Int, Vector{S}}()
    NFcell = length(first(mesh.faces))              # nodes per cell

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

1D-in-2D circle boundary as a CurvilinearMesh with line elements of order `porder`.
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
        gmsh.option.setNumber("Mesh.SecondOrderLinear",  porder >= 2 ? 0 : 1)

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

# ---------------------------
# Small utilities
# ---------------------------
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
sym(A) = 0.5 .* (A .+ transpose(A))
const CHECK_X_EXACT = false  # flip to true to compare bases

# ---------------------------
# Main sweep p = 1:5
# ---------------------------
R = 1.0
N = 64
hs = 2π*R/N

@testset "Curvilinear circle p=1:5 — assembly, SS rules, refinement" begin
    for p in 1:5
        @testset "order p=$p" begin
            # Mesh & bases
            m  = circle_curvilinear(R, p; h=hs)
            X  = lagrangecxd0(m)
            Xe = lagrangec0_gmsh_exact(m)

            # Optional basis equivalence (disable if slow)
            if CHECK_X_EXACT
                Schk = Helmholtz2D.singlelayer(wavenumber=1.0)
                Zx   = sym(assemble(Schk, X,  X;
                                    threading=BEAST.Threading{:single},
                                    quadstrat=BEAST.DoubleNumSauterQstrat(10,10,10,10,10,10)))
                Ze   = sym(assemble(Schk, Xe, Xe;
                                    threading=BEAST.Threading{:single},
                                    quadstrat=BEAST.DoubleNumSauterQstrat(10,10,10,10,10,10)))
                @test norm(Zx - Ze) / max(norm(Ze), eps()) < 1e-10
            end

            # Assembly & symmetry
            S = Helmholtz2D.singlelayer(wavenumber=1.0)
            Zq  = assemble(S, X, X;
                           threading = BEAST.Threading{:single},
                           quadstrat = BEAST.DoubleNumQStrat(15,14))
            @test norm(Zq - transpose(Zq)) / max(norm(Zq), eps()) < 1e-12

            Zss = assemble(S, X, X;
                           threading = BEAST.Threading{:single},
                           quadstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30))
            @test norm(Zss - transpose(Zss)) / max(norm(Zss), eps()) < 1e-12

            δ = norm(sym(Zq) - sym(Zss)) / max(norm(sym(Zss),), eps())
            @test δ < 5e-3

            # Contact classification & rule selection
            is = first(cells(m))
            it = is == last(cells(m)) ? first(cells(m)) : is + 1
            τ  = chart(m, is)
            σ  = chart(m, it)
            k  = mod1(is + max(5, length(cells(m)) ÷ 3), length(cells(m)))
            ρ  = chart(m, k)

            @test edgehits_local(τ, τ) == 2
            @test edgehits_local(τ, σ) == 1
            @test edgehits_local(τ, ρ) == 0

            g  = BEAST.refspace(X)
            qs = BEAST.DoubleNumSauterQstrat(10,10, 10,10, 10,10)
            qd = BEAST.quaddata(S, g, g, [τ], [σ], qs)

            @test BEAST.quadrule(S, g, g, 1, τ, 1, τ, qd, qs) isa SSQ1D.CommonEdge
            @test BEAST.quadrule(S, g, g, 1, τ, 1, σ, qd, qs) isa SSQ1D.CommonVertex
            let r = BEAST.quadrule(S, g, g, 1, τ, 1, ρ, qd, qs)
                @test !(r isa SSQ1D.CommonEdge) && !(r isa SSQ1D.CommonVertex)
            end

            # Quadrature refinement (convergence towards higher-order ref)
            ref = assemble(S, X, X;
                           threading = BEAST.Threading{:single},
                           quadstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30))
            for q in (10, 15, 20)
                Z = assemble(S, X, X;
                             threading = BEAST.Threading{:single},
                             quadstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,q,q))
                δq = norm(Z - ref) / max(norm(ref), eps())
                @test δq < 5e-3
            end
        end
    end
end
