struct Integrand{Op,LSt,LSb,Elt,Elb}
    operator::Op
    local_test_space::LSt
    local_trial_space::LSb
    test_chart::Elt
    trial_chart::Elb
end


function (igd::Integrand)(u,v)
    
    x = neighborhood(igd.test_chart,u)
    y = neighborhood(igd.trial_chart,v)
    
    f = igd.local_test_space(x)
    g = igd.local_trial_space(y)

    return jacobian(x) * jacobian(y) * igd(x,y,f,g)
end

getvalue(a::SVector{N}) where {N} = SVector{N}(getvalue(a.data))
getvalue(a::NTuple{1}) = (a[1].value,)
getvalue(a::NTuple{N}) where {N} = tuple(a[1].value, getvalue(Base.tail(a))...)

getdivergence(a::SVector{N}) where {N} = SVector{N}(getdivergence(a.data))
getdivergence(a::NTuple{1}) = (a[1].divergence,)
getdivergence(a::NTuple{N}) where {N} = tuple(a[1].divergence, getdivergence(Base.tail(a))...)

function _krondot_gen(a::Type{U}, b::Type{V}) where {U<:SVector{N}, V<:SVector{M}} where {M,N}
    ex = :(SMatrix{N,M}(()))
    for m in 1:M
        for n in 1:N
            push!(ex.args[2].args, :(dot(a[$n], b[$m])))
        end
    end
    return ex
end
@generated function _krondot(a::SVector{N}, b::SVector{M}) where {M,N}
    ex = _krondot_gen(a,b)
    return ex
end

function _integrands_gen(::Type{U}, ::Type{V}) where {U<:SVector{N}, V<:SVector{M}} where {M,N}
    ex = :(SMatrix{N,M}(()))
    for m in 1:M
        for n in 1:N
            # push!(ex.args[2].args, :(dot(a[$n], b[$m])))
            push!(ex.args[2].args, :(f(a[$n], b[$m])))
        end
    end
    return ex
end
@generated function _integrands(f, a::SVector{N}, b::SVector{M}) where {M,N}
    ex = _integrands_gen(a,b)
    # println(ex)
    return ex
end




function _integrands_leg_gen(f::Type{U}, g::Type{V}) where {U<:SVector{N}, V<:SVector{M}} where {M,N}
    ex = :(SMatrix{N,M}(()))
    for m in 1:M
        for n in 1:N
            push!(ex.args[2].args, :(integrand(op, kervals, f[$n], x, g[$m], y)))
        end
    end
    return ex
end
@generated function _integrands_leg(op, kervals, f::SVector{N}, x, g::SVector{M}, y) where {M,N}
    _integrands_leg_gen(f, g)
end


# Support for legacy kernels
function (igd::Integrand)(x,y,f,g)

    op = igd.operator
    kervals = kernelvals(op, x, y)
    _integrands_leg(op, kervals, f, x, g, y)

end


function momintegrals!(op::Operator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_chart, trial_chart, out, rule::SauterSchwabStrategy)

    I, J, _, _ = SauterSchwabQuadrature.reorder(
        test_chart.vertices,
        trial_chart.vertices, rule)

    # permute_vertices reparametrizes the simplex without affecting the normal
    test_chart = CompScienceMeshes.permute_vertices(test_chart, I)
    trial_chart = CompScienceMeshes.permute_vertices(trial_chart, J)

    igd = Integrand(op, test_local_space, trial_local_space, test_chart, trial_chart)
    G = SauterSchwabQuadrature.sauterschwab_parameterized(igd, rule)

    QTest = dof_perm_matrix(test_local_space, I)
    QTrial = dof_perm_matrix(trial_local_space, J)
    out_temp = zeros(eltype(out), numfunctions(test_local_space),numfunctions(trial_local_space))
    out_temp = QTest*G*QTrial'
    out[1:numfunctions(test_local_space),1:numfunctions(trial_local_space)] = out_temp

    nothing
end


function momintegrals_test_refines_trial!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    quadrule, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    momintegrals!(op, test_local_space, trial_local_space,
        test_chart, trial_chart, out, quadrule)
end

# const MWOperator3D = Union{MWSingleLayer3D, MWDoubleLayer3D}
function momintegrals_test_refines_trial!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    qr::SauterSchwabStrategy, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    test_mesh = geometry(test_functions)
    trial_mesh = geometry(trial_functions)

    parent_mesh = CompScienceMeshes.parent(test_mesh)
    trial_charts = [chart(test_mesh, p) for p in CompScienceMeshes.children(parent_mesh, trial_cell)]

    qd = quaddata(op, test_local_space, trial_local_space,
        [test_chart], trial_charts, quadstrat)

    for (q,chart) in enumerate(trial_charts)
        qr = quadrule(op, test_local_space, trial_local_space,
            1, test_chart, q ,chart, qd, quadstrat)

        Q = restrict(trial_local_space, trial_chart, chart)
        zlocal = zero(out)

        momintegrals!(op, test_local_space, trial_local_space,
            test_chart, chart, zlocal, qr)

        for j in 1:numfunctions(trial_local_space)
            for i in 1:numfunctions(test_local_space)
                for k in 1:size(Q, 2)
                    out[i,j] += zlocal[i,k] * Q[j,k]
end end end end end



function momintegrals_trial_refines_test!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    quadrule, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    momintegrals!(op, test_local_space, trial_local_space,
        test_chart, trial_chart, out, quadrule)
end


function momintegrals_trial_refines_test!(out, op,
    test_functions, test_cell, test_chart,
    trial_functions, trial_cell, trial_chart,
    qr::SauterSchwabStrategy, quadstrat)

    test_local_space = refspace(test_functions)
    trial_local_space = refspace(trial_functions)

    test_mesh = geometry(test_functions)
    trial_mesh = geometry(trial_functions)

    parent_mesh = CompScienceMeshes.parent(trial_mesh)
    test_charts = [chart(trial_mesh, p) for p in CompScienceMeshes.children(parent_mesh, test_cell)] 

    qd = quaddata(op, test_local_space, trial_local_space,
        test_charts, [trial_chart], quadstrat)

    for (p,chart) in enumerate(test_charts)
        qr = quadrule(op, test_local_space, trial_local_space,
            p, chart, 1, trial_chart, qd, quadstrat)

        Q = restrict(test_local_space, test_chart, chart)
        zlocal = zero(out)
        momintegrals!(op, test_local_space, trial_local_space,
            chart, trial_chart, zlocal, qr)

        for j in 1:numfunctions(trial_local_space)
            for i in 1:numfunctions(test_local_space)
                for k in 1:size(Q, 2)
                    out[i,j] += Q[i,k] * zlocal[k,j]
        end end end
    end
end
