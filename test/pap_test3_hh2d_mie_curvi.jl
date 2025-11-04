using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions
using Gmsh


import CompScienceMeshes: measure, nodes, Simplex
import LinearAlgebra: norm
import CompScienceMeshes: indices, CurvilinearMesh, chart

const Shape = BEAST.Shape
const LagrangeBasis = BEAST.LagrangeBasis


# measure for straight 1D edges embedded in R^U
function measure(s::Simplex{U,1,<:Any,<:Any,T}) where {U,T}
    p0, p1 = nodes(s)
    return norm(p1 - p0)
end

#=
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
=#

function lagrangec0_curvilinear(mesh::CurvilinearMesh; dirichlet::Bool=true)
    VT = eltype(mesh.vertices)
    T  = eltype(VT)
    S  = BEAST.Shape{T}

    faces = mesh.faces
    @assert !isempty(faces) "Mesh has no faces"

    Nloc = length(first(faces))
    @assert all(length(f) == Nloc for f in faces) "Mixed-order faces not supported"
    p = Nloc - 1
    @assert p ≥ 1 "Expecting p ≥ 1 (got Nloc = $Nloc)"

    _end_interior_indices = function (ζnodes; atol = 10*eps(Float64))
        ζ = collect(ζnodes)
        is0(i) = isapprox(ζ[i], 0.0; atol=atol, rtol=0.0)
        is1(i) = isapprox(ζ[i], 1.0; atol=atol, rtol=0.0)
        iL = findfirst(is0, eachindex(ζ)); iL === nothing && (iL = 1)
        iR = findfirst(is1, eachindex(ζ)); iR === nothing && (iR = length(ζ))
        interior = [i for i in eachindex(ζ) if i != iL && i != iR]
        return (iL, iR, interior)
    end

    endpoint_dofs = Dict{Int, Vector{S}}()
    pos = Vector{VT}()
    fns = Vector{Vector{S}}()

    nv    = length(mesh.vertices)
    deg   = zeros(Int, nv)
    is_ep = falses(nv)

    @inbounds for (c, conn) in pairs(faces)
        ch = chart(mesh, conn)
        iL, iR, interior = _end_interior_indices(ch.ζnodes)

        for r in interior
            push!(pos, mesh.vertices[conn[r]])
            push!(fns, [ S(c, r, one(T)) ])
        end

        gL, gR = conn[iL], conn[iR]
        push!(get!(endpoint_dofs, gL, Vector{S}()), S(c, iL, one(T)))
        push!(get!(endpoint_dofs, gR, Vector{S}()), S(c, iR, one(T)))
        deg[gL] += 1; deg[gR] += 1
        is_ep[gL] = true; is_ep[gR] = true
    end

    is_closed_loop = all(i -> (!is_ep[i]) || deg[i] == 2, 1:nv)
    kept_endpoints =
        dirichlet ? collect(keys(endpoint_dofs)) :
        (is_closed_loop ? Int[] : [g for g in keys(endpoint_dofs) if deg[g] == 1])

    sort!(kept_endpoints)
    for g in kept_endpoints
        push!(pos, mesh.vertices[g])
        push!(fns, endpoint_dofs[g])
    end

    NF = Nloc
    return BEAST.LagrangeBasis{p, 0, NF}(mesh, fns, pos)  # <— qualify with BEAST.
end








function lagrangec0d1_curvilinear_exact_linear(mesh::CurvilinearMesh)
    T = coordtype(mesh)
    S = BEAST.Shape{T}

    # Corner-only vertex DOFs (ignore mid nodes on curved edges)
    dofs = Dict{Int, Vector{S}}()
    for (c, cell) in pairs(mesh.faces)
        @assert length(cell) ≥ 2
        g1, g2 = cell[1], cell[2]   # corners (s,e)
        push!(get!(dofs, g1, Vector{S}()), S(c, 1, T(1)))
        push!(get!(dofs, g2, Vector{S}()), S(c, 2, T(1)))
    end

    gids = sort!(collect(keys(dofs)))
    fns  = [dofs[g] for g in gids]
    pos  = [mesh.vertices[g] for g in gids]

    return BEAST.LagrangeBasis{1, 0, 2}(mesh, fns, pos)
end

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



# This unit tests partially repeats hh2d_nearfield, but it
# tests the excitation part more in detail.
# Overall, it is rather an integration test and maybe should
# be part of an example file.
#let
    ε0 = 8.854187821e-12
    μ0 = 4π*1e-7
    c0 = 1/sqrt(ε0*μ0)
    η0 = sqrt(μ0/ε0)

    f = 1e9 # 1 GHz
    ω = 2π*f
    λ = c0/f
    k = 2π/λ
    h = λ/10
    a = 1.0 # radius of the scatterer
    p_order = 2;

    #circle = CompScienceMeshes.meshcircle(a, h)
    #X0 = lagrangecxd0(circle)
    #X0 = lagrangec0d1(circle)

    circle = circle_curvilinear(a, p_order;h)
    #X0 = lagrangec0_gmsh_exact(circle)
    
    X0 = lagrangec0d1_curvilinear_exact_linear(circle)
    #X0 = lagrangec0_curvilinear(circle)

    # Computing the fields on a circle of radius rc
    #rc = 500.0*a
    #pts = meshcircle(rc, 0.3 * rc).vertices

    rc = 500.0*a
    hps = 0.3 * rc
    pts = circle_curvilinear(rc, p_order;h=hps).vertices

    ## Analytical solutions

    dbesselj(n,x) = (besselj(n-1, x)  - besselj(n+1, x)) / 2
    dhankelh2(n,x) = (hankelh2(n-1, x)  - hankelh2(n+1, x)) / 2

    cart2polar(x,y) = SVector(sqrt(x^2 + y^2), atan(y, x))




    ## Analytical Solutions

        # TM planewave solution (based on Sec 6.4.2, Jin, Theory and Computation of Electromagnetic Fields])
        function TM_pec_planewave(E0, k, a, ρ, φ)
            n = 0

            # Jin (6.4.11)
            a_n(n) = -(1.0im)^(-n)*(besselj(n, k*a) / hankelh2(n, k*a))
            valEz(n) =  a_n(n) * hankelh2(n, k*ρ)*exp(im*n*φ)
            valHφ(n) = 1/(im * ω * μ0) * k * a_n(n) * dhankelh2(n,k*ρ) * exp(im*n*φ)
            valHρ(n) = -1/(im * ω * μ0)  * (im * n / ρ) * a_n(n) * hankelh2(n, k*ρ) * exp(im*n*φ)

            retEz = valEz(0)
            retHφ = valHφ(0)
            retHρ = valHρ(0)

            while true
                n += 1
                bufEz = valEz(n) + valEz(-n)
                if abs(bufEz) >= abs(retEz) * eps(eltype(real(retEz))) * 1e3
                    retEz += bufEz
                else
                    break
                end
            end

            n=0
            while true
                n += 1
                bufHφ = valHφ(n) + valHφ(-n)
                if abs(bufHφ) >= abs(retHφ) * eps(eltype(real(retHφ))) * 1e3
                    retHφ += bufHφ
                else
                    break
                end
            end

            n=0
            while true
                n += 1
                bufHρ = valHρ(n) + valHρ(-n)
                if abs(bufHρ) >= abs(retHρ) * eps(eltype(real(retHρ))) * 1e3
                    retHρ += bufHρ
                else
                    break
                end
            end

            return E0 * retEz, E0*(retHρ*SVector(cos(φ), sin(φ)) + retHφ*SVector(-sin(φ), cos(φ)))
        end

        function TM_pec_planewave_E(E0, k, a, pts)
            return [TM_pec_planewave(E0, k, a, cart2polar(p[1], p[2])...)[1] for p in pts]
        end

        function TM_pec_planewave_H(E0, k, a, pts)
            return [TM_pec_planewave(E0, k, a, cart2polar(p[1], p[2])...)[2] for p in pts]
        end


        # TE planewave solution
        function TE_pec_planewave(H0, k, a, ρ, φ)
            n = 0

            # Jin (6.4.19)
            b(n) = -(1.0im)^(-n)*dbesselj(n,k*a)/dhankelh2(n,k*a)

            # Hz is not needed, but maybe we use it later for MFIE
            valHz(n) = b(n) * hankelh2(n, k*ρ)*exp(im*n*φ)
            valEφ(n) = -1/(im * ω * ε0) * k * b(n) * dhankelh2(n, k*ρ)*exp(im*n*φ)
            valEρ(n) = 1/(im * ω * ε0)  * (im * n / ρ) * b(n) * hankelh2(n, k*ρ)*exp(im*n*φ)

            retHz = valHz(0)
            retEφ = valEφ(0)
            retEρ = valEρ(0)

            while true
                n += 1
                bufHz = valHz(n) + valHz(-n)
                if abs(bufHz) >= abs(retHz) * eps(eltype(real(retHz))) * 1e3
                    retHz += bufHz
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufEφ = valEφ(n) + valEφ(-n)
                if abs(bufEφ) >= abs(retEφ) * eps(eltype(real(retEφ))) * 1e3
                    retEφ += bufEφ
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufEρ = valEρ(n) + valEρ(-n)
                if abs(bufEρ) >= abs(retEρ) * eps(eltype(real(retEρ))) * 1e3
                    retEρ += bufEρ
                else
                    break
                end
            end

            return H0*retHz, H0*(retEρ*SVector(cos(φ), sin(φ)) + retEφ*SVector(-sin(φ), cos(φ)))
        end

        function TE_pec_planewave_H(H0, k, a, pts)
            return [TE_pec_planewave(H0, k, a, cart2polar(p[1],p[2])...)[1] for p in pts]
        end

        function TE_pec_planewave_E(H0, k, a, pts)
            return [TE_pec_planewave(H0, k, a, cart2polar(p[1],p[2])...)[2] for p in pts]
        end


        #=
        This follows Section 6.5.3, Jin, Theory and Computation of EM fields. Particulary, Equation 6.5.38 is
        used to obtain the scattered fields.
        =#
        function TM_pec_line_curr(I, k, a, ρ, φ, ρp, φp)

            A_d(n) = besselj(n,k*a)/hankelh2(n,k*a)

            valEz(n) = A_d(n) * hankelh2(n,k*ρ) * hankelh2(n,k*ρp) * exp(im*n*(φ-φp))
            valHρ(n) = -1/(im * ω * μ0)  * (im * n / ρ) * A_d(n) * hankelh2(n,k*ρp) * hankelh2(n,k*ρ) * exp(im*n*(φ-φp))
            valHφ(n) = 1/(im * ω * μ0) * k * A_d(n) * hankelh2(n,k*ρp) * dhankelh2(n,k*ρ) * exp(im*n*(φ-φp))

            retEz = valEz(0)
            retHρ = valHρ(0)
            retHφ = valHφ(0)

            n = 0
            while true
                n += 1
                bufEz = valEz(n) + valEz(-n)
                if abs(bufEz) >= abs(retEz) * eps(eltype(real(retEz))) * 1e3
                    retEz += bufEz
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufHρ = valHρ(n) + valHρ(-n)
                if abs(bufHρ) >= abs(retHρ) * eps(eltype(real(retHρ))) * 1e3
                    retHρ += bufHρ
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufHφ = valHφ(n) + valHφ(-n)
                if abs(bufHφ) >= abs(retHφ) * eps(eltype(real(retHφ))) * 1e3
                    retHφ += bufHφ
                else
                    break
                end
            end

            Hₜ = retHρ * SVector(cos(φ),sin(φ)) + retHφ * SVector(-sin(φ),cos(φ))

            return η0*k*I/4 * retEz, η0*k*I/4 * Hₜ
        end

        function TM_pec_line_curr_E(I, k, a, pts, spt)
            return [TM_pec_line_curr(I,k,a,cart2polar(p[1],p[2])..., spt[1],spt[2])[1] for p in pts]
        end

        function TM_pec_line_curr_H(I, k, a, pts, spt)
            return [TM_pec_line_curr(I,k,a,cart2polar(p[1],p[2])..., spt[1],spt[2])[2] for p in pts]
        end


        #=
        Here, the current source is also infinitely long along the z-direction, but the current is directed in
        the ϕ-direction. The expressions for the scattered field are derived stating that the Hz field satisfies
        the Helmholtz equation and thus can be expanded in terms of harmonics.
        =#

        function An_TE(a::F, n::I; Z = ComplexF64(0)) where {F,I}
            x = k*a
            # TODO: extend to include the impedance BC part
            An_te = dhankelh2(n,x)/dbesselj(n,x)
            return ComplexF64(-An_te)
        end

        function TE_pec_line_curr(I::F, k::F, a::F,
            ρ::F, φ::F, ρp::F, φp::F; 
            Z = ComplexF64(0)) where F

            valHz(n) = dhankelh2(n,k*ρp)/An_TE(a,n) * hankelh2(n,k*ρ) * exp(im*n*(φ-φp))
            valEρ(n) = -1 / (im * ω * ε0)  * (im * n / ρ) * dhankelh2(n,k*ρp) / An_TE(a,n) * hankelh2(n,k*ρ) * exp(im*n*(φ-φp))
            valEφ(n) =  1 / (im * ω * ε0)  * k *  dhankelh2(n,k*ρp) / An_TE(a,n) * dhankelh2(n,k*ρ) * exp(im*n*(φ-φp))

            retHz = valHz(0)
            retEρ = valEρ(0)
            retEφ = valEφ(0)

            n = 0
            while true
                n += 1
                bufHz = valHz(n) + valHz(-n)

                if abs(bufHz) >= abs(retHz) * eps(eltype(real(retHz))) * 1e3
                    retHz += bufHz
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufEρ = valEρ(n) + valEρ(-n)

                if abs(bufEρ) >= abs(retEρ) * eps(eltype(real(retEρ))) * 1e3
                    retEρ += bufEρ
                else
                    break
                end
            end

            n = 0
            while true
                n += 1
                bufEφ = valEφ(n) + valEφ(-n)

                if abs(bufEφ) >= abs(retEφ) * eps(eltype(real(retEφ))) * 1e3
                    retEφ += bufEφ
                else
                    break
                end
            end

            prefac = k^2/(μ0)
            return  prefac * -I*k/(4*im) * retHz, I*k^3/(4*im*μ0) * (retEρ * SVector(cos(φ),sin(φ)) + retEφ * SVector(-sin(φ),cos(φ)))
        end

        function TE_pec_line_curr_H(I, k, a, pts, spt)
            return [TE_pec_line_curr(I,k,a,cart2polar(p[1],p[2])..., spt[1],spt[2])[1] for p in pts]
        end

        function TE_pec_line_curr_E(I, k, a, pts, spt)
            return [TE_pec_line_curr(I,k,a,cart2polar(p[1],p[2])..., spt[1],spt[2])[2] for p in pts]
        end


    ## TM Scattering

        # TM-EFIE system matrix
        𝒮 = Helmholtz2D.singlelayer(; alpha=im*k*η0, wavenumber=k)
        M_TMEFIE = assemble(𝒮, X0, X0)

        # TM-MFIE system matrix
        𝒟ᵀ = Helmholtz2D.doublelayer_transposed(; wavenumber=k)
        Dᵀ = assemble(𝒟ᵀ, X0, X0)

        I0 = assemble(BEAST.Identity(), X0, X0)
        M_TMMFIE = +0.5*I0 - Dᵀ

        # 1. Excitation: Plane wave (The incoming planewave is along the x-axis)
        E0 = 1.0 # Amplitude
        Ez_pw_inc = Helmholtz2D.planewave(; amplitude=E0, wavenumber=k, direction=SVector(1.0, 0.0))

        ez_pw_inc = assemble(DirichletTrace(Ez_pw_inc), X0)
        j_TMEFIE_pw =  M_TMEFIE \ ez_pw_inc

        # TM-EFIE: We compute the scattered Ez component (scalar)
        Ez_pw_sca_num = -potential(HH2DSingleLayerNear(𝒮), pts, j_TMEFIE_pw, X0; type=ComplexF64)
        Ez_pw_sca_ana = TM_pec_planewave_E(E0, k, a, pts)

        @test norm(Ez_pw_sca_ana - Ez_pw_sca_num) ./ norm(Ez_pw_sca_ana) < 0.002

        Ht_pw_sca_num = -potential(HH2DDoubleLayerTransposedNear(𝒟ᵀ), pts, j_TMEFIE_pw, X0; type=SVector{2,ComplexF64})
        Ht_pw_sca_ana = TM_pec_planewave_H(E0, k, a, pts)

        @test norm(Ref(ẑ) .× Ht_pw_sca_num - Ht_pw_sca_ana) ./ norm(Ht_pw_sca_ana) < 0.002

        # TM-MFIE
        Ht_pw_inc =  - 1.0 / (im * ω * μ0) * curl(Ez_pw_inc)
        ht_pw_inc = -assemble(TangentTrace(Ht_pw_inc), X0)

        j_TMMFIE_pw =  M_TMMFIE \ -ht_pw_inc

        @test norm(j_TMMFIE_pw - j_TMEFIE_pw) / norm(j_TMEFIE_pw) < 0.05
        

        #j_TMMFIE_pw__ =  M_TMEFIE \ ht_pw_inc
        #
        #norm(j_TMMFIE_pw__ - j_TMEFIE_pw) / norm(j_TMEFIE_pw)

        # 2. Excitation: Infinite line current
        I = 1.0 # 1 A
        ρp = 2.0 # Source radial pos
        φp = 0.0 # Source angular pos
        xp = ρp * cos(φp)
        yp = ρp * sin(φp)

        # i) Ez Field
        # Choosing the amplitude of Einc such that [Eq 6.5.11, Jin] is satisfied
        Ez_lc_Einc = Helmholtz2D.monopole(;position = SVector(xp,yp), wavenumber=k, amplitude = -η0*k*I*im)
        ez_lc_Einc = assemble(DirichletTrace(Ez_lc_Einc),X0)
        j_TMEFIE_lc = M_TMEFIE \ ez_lc_Einc

        Ez_lc_sca_num = -potential(HH2DSingleLayerNear(𝒮),pts , j_TMEFIE_lc, X0;type=ComplexF64)
        Ez_lc_sca_ana = TM_pec_line_curr_E(I, k, a, pts, SVector(ρp,φp))

        @test norm(Ez_lc_sca_num - Ez_lc_sca_ana) / norm(Ez_lc_sca_ana) < 0.003

        # ii) Hx, Hy Fields
        Ht_lc_sca_num = -potential(HH2DDoubleLayerTransposedNear(𝒟ᵀ), pts, j_TMEFIE_lc, X0; type=SVector{2, ComplexF64})
        Ht_lc_sca_ana = TM_pec_line_curr_H(I, k, a, pts, SVector(ρp,φp))
        @test norm(Ref(ẑ) .× Ht_lc_sca_num - Ht_lc_sca_ana) ./ norm(Ht_lc_sca_ana) < 0.003

        # usage of the TMMFIE current results in a loss of accuracy
        Ht_lc_inc = -1.0 / (im * ω * μ0) * curl(Ez_lc_Einc)
        𝗵t_lc = assemble(TangentTrace(Ht_lc_inc), X0)
        𝗷_TMMFIE_lc = M_TMMFIE \ 𝗵t_lc

        Ht_lc_sca_num = -potential(HH2DDoubleLayerTransposedNear(𝒟ᵀ), pts, 𝗷_TMMFIE_lc, X0; type=SVector{2, ComplexF64})
        norm(Ref(ẑ) .× Ht_lc_sca_num - Ht_lc_sca_ana) ./ norm(Ht_lc_sca_ana)

        # Additional Observation: error increases as the distance between r' and r reduces, i.e., when we make rc smaller.




    ## TE Scattering
        #=
         Incident Field

         Let's first compare the incident fields. To do so, we place the source at ρ = 2.0 and φ = 0.0, and a
         circle of radius r_inc = 1.0 with it's centre at the origin. We compute the incident magnetic fields on the
         circumference of the circle, both analytically and numerically and compare the two, which, naturally, are
         supposed to match.
        =#
        I = 1.0
        ρp = 2.0
        φp = 0.0
        r_inc = 1.0

        pts_inc =  meshcircle(r_inc, 0.6*r_inc).vertices

        function Hz_inc(I,k,ρ,φ,ρp,φp)

            val(n) = dhankelh2(n,k*ρp) * besselj(n,k*ρ) * exp(im*n*φ)
            ret = val(0)

            n = 0
            while true
                n += 1
                buf = val(n) + val(-n)
                if abs(buf) >= abs(ret) * eps(eltype(real(ret))) * 1e3
                    ret += buf
                else
                    break
                end
            end
            prefac = k^2/(μ0)
            return  prefac * -I*k/(4*im) * ret
        end

        function Hz_inc(I::F, k::F, pts::Vector{SVector{2,F}}, spt::SVector{2,F}) where F
            return [Hz_inc(I, k, cart2polar(p[1],p[2])..., spt[1], spt[2]) for p in pts]
        end

        xp = ρp * cos(φp)
        yp = ρp * sin(φp)

        # amplitude = k^2/mu0 to make it a magnetic monopole
        Hz_lc_inc_num = Helmholtz2D.directedmonopole(;position = SVector(xp, yp),
            direction = SVector(0.0,1.0),
            wavenumber = k,
            amplitude = k^2/μ0)
        Hz_lc_inc_ana = Hz_inc(I, k, pts_inc, SVector(ρp, φp))

        @test norm( Hz_lc_inc_ana - Hz_lc_inc_num.(pts_inc)) /norm(Hz_lc_inc_ana) ≈ 0 atol = 1e-12

        # Scattered Field
        # TE-MFIE
        𝒟 = Helmholtz2D.doublelayer(; wavenumber=k)
        D = assemble(𝒟, X1, X1)
        I1 = assemble(BEAST.Identity(), X1, X1)
        M_TEMFIE = 0.5I1 - D

        # TE-EFIE
        𝒩 = Helmholtz2D.hypersingular(;alpha=im*ω*μ0, beta=1.0/(im*ω*ε0), wavenumber=k)
        M_TEEFIE = assemble(𝒩, X1, X1)

        # 1. Excitation: Planewave
        H0 = 1.0
        Hz_pw_inc = Helmholtz2D.planewave(; amplitude=H0, wavenumber=k, direction=SVector(1.0, 0.0))
        hz_pw_inc = assemble(DirichletTrace(Hz_pw_inc), X1)
        j_TEMFIE_pw = M_TEMFIE \ hz_pw_inc

        Hz_pw_sca_num = potential(HH2DDoubleLayerNear(𝒟), pts, j_TEMFIE_pw, X1; type=ComplexF64)
        Hz_pw_sca_ana = TE_pec_planewave_H(H0, k, a, pts)

        @test norm(Hz_pw_sca_num - Hz_pw_sca_ana) /norm(Hz_pw_sca_ana) < 0.0015

        Et_pw_inc = 1 / (im * ω * ε0) * curl(Hz_pw_inc)
        @test Et_pw_inc.direction == Hz_pw_inc.direction
        @test Et_pw_inc.polarization == SVector(0.0, 1.0)
        @test Et_pw_inc.amplitude ≈ Et_pw_inc.gamma / (im * ω * ε0) atol=1e-15
        et_pw_inc = assemble(TangentTrace(Et_pw_inc), X1)

        j_TEEFIE_pw = M_TEEFIE \ et_pw_inc

        @test norm(j_TEEFIE_pw - j_TEMFIE_pw)/norm(j_TEMFIE_pw) < 0.007

        Et_pw_sca_num = -potential(HH2DHyperSingularNear(𝒩), pts, j_TEEFIE_pw, X1; type=SVector{2, ComplexF64})
        Et_pw_sca_ana = TE_pec_planewave_E(H0, k, a, pts)

        # We compute the scattered Ez component (scalar)
        @test norm(Ref(ẑ) .× Et_pw_sca_num - Et_pw_sca_ana) / norm(Et_pw_sca_ana) < 0.003

        # 2. Excitation: Infinite line current
        # i) Hz field
        hz_lc_inc = -assemble(DirichletTrace(Hz_lc_inc_num),X1)
        j_TEMFIE_lc = M_TEMFIE \ hz_lc_inc

        Hz_lc_sca_ana = TE_pec_line_curr_H(I, k, a, pts, SVector(ρp, φp))

        Hz_lc_sca_num = -potential(HH2DDoubleLayerNear(𝒟), pts, j_TEMFIE_lc, X1; type=ComplexF64)
        @test norm( Hz_lc_sca_ana - Hz_lc_sca_num) /norm(Hz_lc_sca_ana) <= 0.003

        # ii) Ex, Ey fields
        # a. RHS with j_TEMFIE_lc
        Et_lc_sca_num = potential(HH2DHyperSingularNear(𝒩), pts, j_TEMFIE_lc, X1; type=SVector{2, ComplexF64})
        Et_lc_sca_num_curl = Ref(ẑ) .× Et_lc_sca_num

        Et_lc_sca_ana = TE_pec_line_curr_E(I, k, a, pts, SVector(ρp, φp))

        @test norm( Et_lc_sca_num_curl - Et_lc_sca_ana) ./ norm(Et_lc_sca_ana) <= 0.003

        # b. RHS with j_TEEFIE
        Et_lc_inc_num = 1.0/(im*ω*ε0) * curl(Hz_lc_inc_num)
        et_lc_inc = -assemble(TangentTrace(Et_lc_inc_num),X1)
        j_TEEFIE_lc = M_TEEFIE \  et_lc_inc

        @test norm(j_TEEFIE_lc - j_TEMFIE_lc) ./ norm(j_TEMFIE_lc) <= 0.01

        Et_lc_sca_num_2 = potential(HH2DHyperSingularNear(𝒩), pts, j_TEEFIE_lc, X1; type=SVector{2, ComplexF64})
        Et_lc_sca_num_2_curl = Ref(ẑ) .× Et_lc_sca_num_2
        @test norm(Et_lc_sca_num_2_curl - Et_lc_sca_ana) ./ norm(Et_lc_sca_ana) <= 0.003
#end