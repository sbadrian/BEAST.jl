struct DoubleNumSauterTellesQstrat{R,T,S} <: AbstractQuadStrat
    outer_rule::R
    inner_rule::R
    telles_outer_rule::T
    telles_inner_rule::T
    sauter_schwab_common_edge::S
    sauter_schwab_common_vert::S
end

function momintegrals!(op::Operator,
    test_local_space, trial_local_space,
    test_chart, trial_chart,
    out, rule::BEAST.TellesQuadrature1D.TellesStrategy1D)
    #rule::Union{SauterSchwabStrategy,SauterSchwabQuadrature1D.SauterSchwabStrategy1D})

    num_tshapes = numfunctions(test_local_space, domain(test_chart))
    num_bshapes = numfunctions(trial_local_space, domain(trial_chart))

    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    #G = Matrix{scalartype(op)}(undef, num_tshapes, num_bshapes)
    #G = zeros(scalartype(op), num_tshapes, num_bshapes)
    G = BEAST.TellesQuadrature1D.telles_parametrized(igd, rule, num_tshapes, num_bshapes)

    for j in 1:num_bshapes
        for i in 1:num_tshapes
            out[i, j] += G[i, j]
        end
    end

    nothing
end
