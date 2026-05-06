module TellesQuadrature1D

# -------- exportet parts
# types
export TellesStrategy1D
# functions
export telles_parametrized
# -------- included files
#include("doubletellesint.jl")

# -------- imports
import ..scalartype

# -------- used packages
using FastGaussQuadrature
using CompScienceMeshes
using LinearAlgebra
using StaticArrays


struct TellesStrategy1D{A}
    qpso::A # GL quadrature for the outer integral
    qpsi::A # GL quadrature for the inner integral
end



const R_BAR_THRESHOLDS = [0.001, 0.005, 0.01, 0.03, 0.075, 0.15, 0.25, 0.4, 0.6, 0.8, 1.25, 2.0, 3.0, 4.0, 4.5]
const R_BAR_VALUES = [0.051, 0.086, 0.116, 0.197, 0.316, 0.446, 0.564, 0.605, 0.731, 0.806, 0.892, 0.946, 0.972, 0.98, 0.99]

function getr_bar(D::Float64)
    if D <= R_BAR_THRESHOLDS[1]
        return R_BAR_VALUES[1]
    elseif D >= R_BAR_THRESHOLDS[end]
        return R_BAR_VALUES[end]
    end

    idx = searchsortedfirst(R_BAR_THRESHOLDS, D)

    # Linear interpolation between idx-1 and idx
    D0 = R_BAR_THRESHOLDS[idx-1]
    D1 = R_BAR_THRESHOLDS[idx]
    R0 = R_BAR_VALUES[idx-1]
    R1 = R_BAR_VALUES[idx]

    return R0 + (D - D0) * (R1 - R0) / (D1 - D0)
end



@inline function telles_deg3(acc, igd, aux, refpoint::SVector, qpsN1)

    # --- TYPE OF ACCUMULATOR ---
    #acc = zero(scalartype(igd.operator))        #type of the result
    Edge = vertices(igd.trial_chart)
    jac = 1 / 2
    qpsN = [(2 * v1 - 1, 2 * w1) for (v1, w1) in qpsN1]
    # --- GEOMETRY ---

    L = volume(igd.trial_chart)
    η_tilde = normalize(igd.trial_chart.tangents[1])
    nullpoint = cartesian(neighborhood(igd.trial_chart, 0.0))
    onepoint = cartesian(neighborhood(igd.trial_chart, 1.0))

    𝒹_tilde = refpoint - nullpoint
    u_star = dot(η_tilde, 𝒹_tilde) / L
    T = nullpoint + η_tilde * L * u_star
    Dist = norm(refpoint - T)
    η_bar = u_star * 2 - 1

    #η_tilde_tilde = normalize((refpoint - T))         #normal vector
    #@show dot(η_tilde_tilde, η_tilde)
    if η_bar > 1
        #η_bar = 1
        D = 2 * norm(refpoint - onepoint) / L
    elseif η_bar < -1
        #η_bar = 0
        D = 2 * norm(refpoint - nullpoint) / L
    else
        D = 2 * Dist / L
    end

    #u_star, D, η_bar, T = getparam2(vertices(igd.trial_chart), refpoint)
    #Dist = norm(refpoint - T)
    #@show T
    #@show dot(T - refpoint, T - Edge[1])
    # --- AXIS CHECK ---
    isonaxis = true
    threshhold = 10^-13
    if Dist > threshhold
        isonaxis = false
    end
    #=
    if abs(η_bar) < 1
        # Recursively split the edge, return type T
        return tellesnocm_integratedeg3(K,
            SVector(nullpoint, nullpoint + u_star * L * η_tilde),
            refpoint, N
        ) + tellesnocm_integratedeg3(K,
            SVector(nullpoint + u_star* L * η_tilde, Edge[1]),
            refpoint, N
        )
    end=#
    # ==============================================================
    # ======================= NOT ON AXIS ===========================
    # ==============================================================
    if !isonaxis

        r_bar = getr_bar(D)

        #η_star = η_bar^2 - 1
        q = (η_bar * (3 - 2 * r_bar) - (2 * η_bar^3) / (1 + 2 * r_bar)) / (2 * (1 + 2 * r_bar)^2) - η_bar / (2 * (1 + 2 * r_bar))
        p = (4 * r_bar * (1 - r_bar) + 3 * (1 - η_bar^2)) / (3 * (1 + 2 * r_bar)^2)

        γ_bar = cbrt(-q + sqrt(q^2 + p^3)) +
                cbrt(-q - sqrt(q^2 + p^3)) +
                η_bar / (1 + 2 * r_bar)

        Q = 1 + 3 * γ_bar^2

        # Polynomial coefficients
        coef = (
            a=(1 - r_bar) / Q,
            b=-3 * (1 - r_bar) * γ_bar / Q,
            c=(r_bar + 3 * γ_bar^2) / Q,
            d=3 * (1 - r_bar) * γ_bar / Q    # (=-b)
        )

        # A let-block keeps mapv / dmapv visible for loop use
        let a = coef.a, b = coef.b, c = coef.c, d = coef.d

            @inline mapv(v1) = a * v1^3 + b * v1^2 + c * v1 + d
            @inline dmapv(v1) = 3 * a * v1^2 + 2 * b * v1 + c

            @inbounds for j in eachindex(qpsN)
                v1, w1 = qpsN[j]

                t = mapv(v1)
                ξ = (t + 1) / 2
                #x = P1 + ξ * Δ

                acc += jac * dmapv(v1) * w1 *
                       aux(ξ)
            end
        end

        return acc
    end

    # ==============================================================
    # ========================== ON AXIS ============================
    # ==============================================================
    η_star = η_bar^2 - 1

    γ_bar = cbrt(η_bar * η_star + abs(η_star)) +
            cbrt(η_bar * η_star - abs(η_star)) +
            η_bar

    denom = 1 + 3 * γ_bar^2

    let γ = γ_bar, denom = denom

        @inline mapv(v1) = ((v1 - γ)^3 + γ * (γ^2 + 3)) / denom
        @inline dmapv(v1) = 3 * (v1 - γ)^2 / denom

        @inbounds for j in eachindex(qpsN)
            v1, w1 = qpsN[j]

            t = mapv(v1)
            ξ = (t + 1) / 2

            #x = P1 + ξ * Δ

            acc += jac * dmapv(v1) * w1 *
                   aux(ξ)
        end
    end

    return acc
end

function telles_parametrized(igd, rule::TellesStrategy1D, num_tshapes, num_bshapes)
    qpsi = rule.qpsi
    qpso = rule.qpso
    G = zeros(scalartype(igd.operator), num_tshapes, num_bshapes)
    for (v1, w1) in qpso
        G += w1 * telles_deg3(zeros(scalartype(igd.operator), num_tshapes, num_bshapes), igd, x -> igd(x, v1), cartesian(igd.test_chart, v1), qpsi)
    end
    return G
end

#=
function sauterschwab_parameterized1D(integrand, strategy::SauterSchwabStrategy1D)
    return sum(w1 * w2 * strategy(integrand, η, ξ) for (η, w1) in strategy.qpsi, (ξ, w2) in strategy.qpso)
end
=#
end
