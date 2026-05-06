using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
#import PlotlyJS       <-- For Plottong only


@testset "Telles-Integration" begin
    ℒ = 10.0
    𝒹 = 2.5 * ℒ / 2
    Stat = Helmholtz2D.singlelayer()
    P1 = SVector(0.0, 0.0)
    P2 = SVector(ℒ, 0.0)
    P3 = SVector(P1[1] + 0.25 * ℒ, P1[2] + 𝒹)
    P4 = SVector(P2[1] + 0.25 * ℒ, P2[2] + 𝒹)
    vertices = [P1, P2, P3, P4]

    el1 = CompScienceMeshes.SimplexGraph{2}(index(1,2))
    el2 = CompScienceMeshes.SimplexGraph{2}(index(3,4))
    faces = [el1, el2]
    tstmesh = Mesh(vertices, faces)
    X01 = lagrangecxd0(tstmesh)
    function referencestat(l::Float64, y::SVector)
        return -(1 / (4 * π)) * (((1 * l - y[1]) * log((1 * l - y[1])^2 + y[2]^2) - 2 * (1 * l - y[1]) + 2 * y[2] * atan((1 * l - y[1]) / y[2])) - ((0 * l - y[1]) * log((0 * l - y[1])^2 + y[2]^2) - 2 * (0 * l - y[1]) + 2 * y[2] * atan((0 * l - y[1]) / y[2])))
    end
    function referencestatE(Edge::SVector{2,<:SVector}, Edge2::SVector{2,<:SVector}, N::Int)
        qpsN = BEAST._legendre(N, -1.0, 1.0)
        J = (norm(Edge2[2] - Edge2[1]) / 2)
        Δx = (Edge2[2][1] - Edge2[1][1])
        Δy = (Edge2[2][2] - Edge2[1][2])
        L = norm(Edge[2] - Edge[1])
        return sum(w1 * J * referencestat(L, SVector(Edge2[1][1] + ((v1 + 1) / 2) * Δx, Edge2[1][2] + ((v1 + 1) / 2) * Δy)) for (v1, w1) in qpsN)
    end
    reforder = 10000
    relerT = Float64[]
    ref = referencestatE(SVector(P1, P2), SVector(P3, P4), reforder)
    quads = [10, 30]
    for i in quads
        quadstrat = BEAST.DoubleNumSauterTellesQstrat(i, i, i, i, i, i)
        tel = assemble(Stat, X01, X01, quadstrat=quadstrat)[1, 2]
        push!(relerT, abs(tel - ref) / abs(ref))
    end


    @test relerT[1] < 10^-3
    @test relerT[2] < 10^-7

    function returntestmesh(𝒹::Float64, α::Float64, ℒ::Float64)
        α = α * π / 180
        P1 = SVector(-ℒ / 2, 0.0)
        P6 = SVector(ℒ / 2, 0.0)
        P3 = SVector(P1[1], P1[2] + 𝒹)
        P4 = SVector(P6[1], P6[2] + 𝒹)
        P2 = SVector(P1[1] - (𝒹 * 0.5 / tan(α / 2)), P1[2] + 𝒹 / 2)
        P5 = SVector(P6[1] + (𝒹 * 0.5 / tan(α / 2)), P6[2] + 𝒹 / 2)
        vertices = [P1, P2, P3, P4, P5, P6]
        el1 = CompScienceMeshes.SimplexGraph{2}(index(1,2))
        el2 = CompScienceMeshes.SimplexGraph{2}(index(2,3))
        el3 = CompScienceMeshes.SimplexGraph{2}(index(3,4))
        el4 = CompScienceMeshes.SimplexGraph{2}(index(4,5))
        el5 = CompScienceMeshes.SimplexGraph{2}(index(5,6))
        el6 = CompScienceMeshes.SimplexGraph{2}(index(6,1))
        faces = [el1, el2, el3, el4, el5, el6]
        M = Mesh(vertices, faces)
        return M
    end
    mesh = returntestmesh(1 / 20, 10.0, 1.0)
    #=========================For Plotting only=========================##=
    function make_3d(M)
        v = M.vertices
        f = M.faces
        new = [SVector(v1[1], v1[2], 0.0) for v1 in v]
        M3d = Mesh(new, f)
        return M3d
    end
    mesh3d = make_3d(mesh)
    plt = PlotlyJS.plot(CompScienceMeshes.wireframe(mesh3d))
    PlotlyJS.relayout!(plt, scene=PlotlyJS.attr(
    aspectmode="data"   # keeps x, y, z scaling equal to data units
    ))
    PlotlyJS.display(plt)
    =##============================End Plots==============================#
    X1 = lagrangec0d1(mesh)
    quads = [10, 30]
    relD = Float64[]
    for i in quads
        quadstrat = BEAST.DoubleNumSauterTellesQstrat(i, i, i, i, i, i)
        default = BEAST.DoubleNumSauterQstrat(i, i, 0, 4, i, i)
        B1 = assemble(Stat, X1, X1, quadstrat=quadstrat)
        B2 = assemble(Stat, X1, X1, quadstrat=default)
        push!(relD, norm(B2 - B1) / norm(B2))
    end
    @test relD[1] < 10^-2
    @test relD[2] < 10^-4

end
