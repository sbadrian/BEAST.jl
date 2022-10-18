using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test

@testset "Helmholtz potential operators" begin
    r = 10.0
    λ = 20 * r
    k = 2 * π / λ

    sphere = meshsphere(r, 0.2 * r)
    X0 = lagrangecxd0(sphere)
    X1 = lagrangec0d1(sphere)

    S = Helmholtz3D.singlelayer(; gamma=im * k)
    D = Helmholtz3D.doublelayer(; gamma=im * k)
    Dt = Helmholtz3D.doublelayer_transposed(; gamma=im * k)
    N = Helmholtz3D.hypersingular(; gamma=im * k)

    q = 100.0
    ϵ = 1.0

    # Interior problem
    # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1
    pos1 = SVector(r * 1.5, 0.0, 0.0)  # positioning of point charges
    pos2 = SVector(-r * 1.5, 0.0, 0.0)

    # Potential of point charges
    function Φ_inc(x)
        return q / (4 * π * ϵ) * (
            exp(-im * k * norm(x - pos1)) / (norm(x - pos1)) -
            exp(-im * k * norm(x - pos2)) / (norm(x - pos2))
        )
    end

    function ∂nΦ_inc(x)
        return -q / (4 * π * ϵ * r) * (
            (dot(
                x,
                (x - pos1) * exp(-im * k * norm(x - pos1)) / (norm(x - pos1)^2) *
                (im * k + 1 / (norm(x - pos1))),
            )) - (dot(
                x,
                (x - pos2) * exp(-im * k * norm(x - pos2)) / (norm(x - pos2)^2) *
                (im * k + 1 / (norm(x - pos2))),
            ))
        )
    end

    gD0 = assemble(ScalarTrace(Φ_inc), X0)
    gD1 = assemble(ScalarTrace(Φ_inc), X1)
    gN = assemble(ScalarTrace(∂nΦ_inc), X1)

    G = assemble(Identity(), X1, X1)
    o = ones(numfunctions(X1))

    # Interior Dirichlet problem
    M_IDPSL = assemble(S, X0, X0) # Single layer (SL)
    M_IDPDL = (-1 / 2 * assemble(Identity(), X1, X1) + assemble(D, X1, X1)) # Double layer (DL)

    # Interior Neumann problem
    # Neumann derivative from DL potential with deflected nullspace
    M_INPDL = -assemble(N, X1, X1) + G * o * o' * G
    # Neumann derivative from SL potential with deflected nullspace
    M_INPSL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1)) + G * o * o' * G 

    ρ_IDPSL = M_IDPSL \ (-gD0)
    ρ_IDPDL = M_IDPDL \ (gD1)

    ρ_INPDL = M_INPDL \ (gN)
    ρ_INPSL = M_INPSL \ (-gN)

    pts = meshsphere(0.8 * r, 0.8 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated

    pot_IDPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_IDPSL, X0; type=ComplexF64)
    pot_IDPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_IDPDL, X1; type=ComplexF64)

    pot_INPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_INPSL, X1; type=ComplexF64)
    pot_INPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_INPDL, X1; type=ComplexF64)

    # Total field inside should be zero
    err_IDPSL_pot = norm(pot_IDPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_IDPDL_pot = norm(pot_IDPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPSL_pot = norm(pot_INPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))
    err_INPDL_pot = norm(pot_INPDL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

    # Analytical solution, based on ???
    function fieldtheo(x)
        return q / (4 * π * ϵ) * (
            (
                (x - pos1) * exp(-im * k * norm(x - pos1)) / (norm(x - pos1)^2) *
                (im * k + 1 / (norm(x - pos1)))
            ) - (
                (x - pos2) * exp(-im * k * norm(x - pos2)) / (norm(x - pos2)^2) *
                (im * k + 1 / (norm(x - pos2)))
            )
        )
    end

    field_IDPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_IDPSL, X0)
    field_IDPDL = potential(HH3DHyperSingularNear(im * k), pts, ρ_IDPDL, X1)
    field_INPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_INPSL, X1)
    field_INPDL = potential(HH3DHyperSingularNear(im * k), pts, ρ_INPDL, X1)

    err_IDPSL_field = norm(field_IDPSL + fieldtheo.(pts)) / norm(fieldtheo.(pts))
    err_IDPDL_field = norm(field_IDPDL + fieldtheo.(pts)) / norm(fieldtheo.(pts))
    err_INPSL_field = norm(field_INPSL + fieldtheo.(pts)) / norm(fieldtheo.(pts))
    err_INPDL_field = norm(field_INPDL + fieldtheo.(pts)) / norm(fieldtheo.(pts))

    # Exterior problem 
    # formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.2

    pos1 = SVector(r * 0.5, 0.0, 0.0)
    pos2 = SVector(-r * 0.5, 0.0, 0.0)

    gD0 = assemble(ScalarTrace(Φ_inc), X0)
    gD1 = assemble(ScalarTrace(Φ_inc), X1)
    gN = assemble(ScalarTrace(∂nΦ_inc), X1)

    G = assemble(Identity(), X1, X1)
    o = ones(numfunctions(X1))

    M_EDPSL = assemble(S, X0, X0)
    M_EDPDL = (1 / 2 * assemble(Identity(), X1, X1) + assemble(D, X1, X1))

    M_ENPDL = -assemble(N, X1, X1) + G * o * o' * G
    M_ENPSL = -1 / 2 * assemble(Identity(), X1, X1) + assemble(Dt, X1, X1) + G * o * o' * G

    ρ_EDPSL = M_EDPSL \ (-gD0)
    ρ_EDPDL = M_EDPDL \ (gD1)

    ρ_ENPDL = M_ENPDL \ gN
    ρ_ENPSL = M_ENPSL \ (-gN)

    testsphere = meshsphere(1.2 * r, 1.2 * 0.6 * r)
    pts = testsphere.vertices[norm.(testsphere.vertices) .> r]

    pot_EDPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_EDPSL, X0; type=ComplexF64)
    pot_EDPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_EDPDL, X1; type=ComplexF64)
    pot_ENPDL = potential(HH3DDoubleLayerNear(im * k), pts, ρ_ENPDL, X1; type=ComplexF64)
    pot_ENPSL = potential(HH3DSingleLayerNear(im * k), pts, ρ_ENPSL, X1; type=ComplexF64)

    err_EDPSL_pot = norm(pot_EDPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_EDPDL_pot = norm(pot_EDPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPSL_pot = norm(pot_ENPSL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))
    err_ENPDL_pot = norm(pot_ENPDL + Φ_inc.(pts)) ./ norm(Φ_inc.(pts))

    field_EDPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_EDPSL, X0)
    field_EDPDL = potential(HH3DHyperSingularNear(im * k), pts, ρ_EDPDL, X1)
    field_ENPSL = -potential(HH3DDoubleLayerTransposedNear(im * k), pts, ρ_ENPSL, X1)
    field_ENPDL = potential(HH3DHyperSingularNear(im * k), pts, ρ_ENPDL, X1)

    err_EDPSL_field = norm(field_EDPSL + fieldtheo.(pts)) / norm(fieldtheo.(pts))
    err_EDPDL_field = norm(field_EDPDL + fieldtheo.(pts)) / norm(fieldtheo.(pts))
    err_ENPSL_field = norm(field_ENPSL + fieldtheo.(pts)) / norm(fieldtheo.(pts))
    err_ENPDL_field = norm(field_ENPDL + fieldtheo.(pts)) / norm(fieldtheo.(pts))

    # errors of interior problems
    @test err_IDPSL_pot < 0.01
    @test err_IDPDL_pot < 0.01
    @test err_INPSL_pot < 0.02
    @test err_INPDL_pot < 0.01

    @test err_IDPSL_field < 0.01
    @test err_IDPDL_field < 0.02
    @test err_INPSL_field < 0.02
    @test err_INPDL_field < 0.01

    # errors of exterior problems
    @test err_EDPSL_pot < 0.01
    @test err_EDPDL_pot < 0.02
    @test err_ENPSL_pot < 0.01
    @test err_ENPDL_pot < 0.01

    @test err_EDPSL_field < 0.01
    @test err_EDPDL_field < 0.03
    @test err_ENPSL_field < 0.01
    @test err_ENPDL_field < 0.02
end