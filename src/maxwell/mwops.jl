abstract type MaxwellOperator3D <: IntegralOperator end
abstract type MaxwellOperator3DReg <: MaxwellOperator3D end

struct KernelValsMaxwell3D{T,U,P,Q}
    "gamma = im * wavenumber"
    gamma::U
    vect::P
    dist::T
    green::U
    gradgreen::Q
end

const inv_4pi = 1/(4pi)
function kernelvals(biop::MaxwellOperator3D, p, q)

    γ = biop.gamma
    r = cartesian(p) - cartesian(q)
    R = norm(r)
    γR = γ*R

    inv_R = 1/R

    expn = exp(-γR)
    green = expn * inv_R * inv_4pi
    gradgreen = -(γ + inv_R) * green * inv_R * r

    KernelValsMaxwell3D(γ, r, R, green, gradgreen)
end

function kernelvals(kernel::MaxwellOperator3DReg, p, q)

    γ = kernel.gamma
    r = p.cart - q.cart
    R = norm(r)
    γR = γ*R

    Exp = exp(-γ*R)
    green = (Exp - 1 + γR - 0.5*γR^2) / (4pi*R)
    gradgreen = ( - (γR + 1)*Exp + (1 - 0.5*γR^2) ) * (r/R^3) / (4π)

    KernelValsMaxwell3D(γ, r, R, green, gradgreen)
end

struct MWSingleLayer3D{T,U} <: MaxwellOperator3D
  gamma::T
  α::U
  β::U
end

MWSingleLayer3D(gamma)  = MWSingleLayer3D(gamma, -gamma, -1/(gamma))
MWWeaklySingular(gamma) = MWSingleLayer3D(gamma, 1, 0)
MWHyperSingular(gamma)  = MWSingleLayer3D(gamma, 0, 1)



export Maxwell3D, MWSingleLayer3DReg

struct MWSingleLayer3DReg{T,U} <: MaxwellOperator3DReg
    gamma::T
    α::U
    β::U
end

struct MWSingleLayer3DSng{T,U} <: MaxwellOperator3D
    gamma::T
    α::U
    β::U
end


function momintegrals!(rop::MWSingleLayer3DReg, g, f, t, s, z, strat::SingularityExtractionRule)
    
  # compute the regular part
  rstrat = regularpart_quadrule(strat)
  momintegrals!(rop, g, f, t, s, z, rstrat)

end

scalartype(op::MaxwellOperator3D) = typeof(op.gamma)

regularpart(op::MWSingleLayer3D) = MWSingleLayer3DReg(op.gamma, op.α, op.β)
singularpart(op::MWSingleLayer3D) = MWSingleLayer3DSng(op.gamma, op.α, op.β)


function _legendre(n,a,b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b-a)/2
    x .= (x.+1)/2*(b-a).+a
    collect(zip(x,w))
end

defaultquadstrat(op::MaxwellOperator3D, tfs::Space, bfs::Space) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)
defaultquadstrat(op::MaxwellOperator3D, tfs::RefSpace, bfs::RefSpace) = DoubleNumWiltonSauterQStrat(2,3,6,7,5,5,4,3)


function quaddata(op::MaxwellOperator3D,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumWiltonSauterQStrat)

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule_far,qs.outer_rule_near))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule_far,qs.inner_rule_near))
    leg = (
      _legendre(qs.sauter_schwab_common_vert,0,1),
      _legendre(qs.sauter_schwab_common_edge,0,1),
      _legendre(qs.sauter_schwab_common_face,0,1),)


    # High accuracy rules (use them e.g. in LF MFIE scenarios)
    # tqd = quadpoints(test_local_space, test_charts, (8,8))
    # bqd = quadpoints(trial_local_space, trial_charts, (8,9))
    # leg = (_legendre(8,a,b), _legendre(10,a,b), _legendre(5,a,b),)


    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end

# use Union type so this code can be shared between the operator
# and its regular part.
MWSL3DGen = Union{MWSingleLayer3D,MWSingleLayer3DReg}
function integrand(biop::MWSL3DGen, kerneldata, tvals, tgeo, bvals, bgeo)

  gx = tvals[1]
  fy = bvals[1]

  dgx = tvals[2]
  dfy = bvals[2]

  G = kerneldata.green
  γ = kerneldata.gamma

  α = biop.α
  β = biop.β

  t = (α * dot(gx, fy) + β * (dgx*dfy)) * G
end

struct MWDoubleLayer3D{T} <: MaxwellOperator3D
  gamma::T
end

struct MWDoubleLayer3DSng{T} <: MaxwellOperator3D
    gamma::T
end

struct MWDoubleLayer3DReg{T} <: MaxwellOperator3DReg
    gamma::T
end

regularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DReg(op.gamma)
singularpart(op::MWDoubleLayer3D) = MWDoubleLayer3DSng(op.gamma)

const MWDL3DGen = Union{MWDoubleLayer3D,MWDoubleLayer3DReg}
function integrand(biop::MWDL3DGen, kerneldata, tvals, tgeo, bvals, bgeo)
    g = tvals[1]
    f = bvals[1]
    ∇G = kerneldata.gradgreen
    (f × g) ⋅ ∇G
end


# quadrule(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd) = qrss(op, g, f, i, τ, j, σ, qd)

# function qrss(op, g, f, i, τ, j, σ, qd)
function quadrule(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace,  i, τ, j, σ, qd,
      qs::DoubleNumWiltonSauterQStrat)

    hits = 0
    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
    dmin2 = floatmax(eltype(eltype(τ.vertices)))
    for t in τ.vertices
        for s in σ.vertices
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < dtol)
        end
    end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    h2 = volume(σ)
    xtol2 = 0.2 * 0.2
    k2 = abs2(op.gamma)
    max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
        qd.tpoints[2,i],
        DoubleQuadRule(
            qd.tpoints[2,i],
            qd.bpoints[2,j],),)
    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end

function quadrule(op::MaxwellOperator3D, g::BDMRefSpace, f::BDMRefSpace,  i, τ, j, σ, qd,
  qs::DoubleNumWiltonSauterQStrat)

  hits = 0
  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
  dmin2 = floatmax(eltype(eltype(τ.vertices)))
  for t in τ.vertices
      for s in σ.vertices
          d2 = LinearAlgebra.norm_sqr(t-s)
          dmin2 = min(dmin2, d2)
          hits += (d2 < dtol)
      end
  end

  hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
  hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
  hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

  h2 = volume(σ)
  xtol2 = 0.2 * 0.2
  k2 = abs2(op.gamma)
  # max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
  #     qd.tpoints[2,i],
  #     DoubleQuadRule(
  #         qd.tpoints[2,i],
  #         qd.bpoints[2,j],),)
  return DoubleQuadRule(
      qd.tpoints[1,i],
      qd.bpoints[1,j],)
end


function qrib(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd)
  # defines coincidence of points
  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

  # decides on whether to use singularity extraction
  xtol = 0.2

  k = norm(op.gamma)

  hits = 0
  xmin = xtol
  for t in τ.vertices
    for s in σ.vertices
      d = norm(t-s)
      xmin = min(xmin, k*d)
      if d < dtol
        hits +=1
        break
      end
    end
  end

    hits == 3   && return BogaertSelfPatchStrategy(5)
    hits == 2   && return BogaertEdgePatchStrategy(8, 4)
    hits == 1   && return BogaertPointPatchStrategy(2, 3)
    rmin = xmin/k
    xmin < xtol && return WiltonSERule(
      qd.tpoints[1,i],
      DoubleQuadRule(
        qd.tpoints[2,i],
        qd.bpoints[2,j],
      ),
    )
    return DoubleQuadRule(
      qd.tpoints[1,i],
      qd.bpoints[1,j],
    )

end


function qrdf(op::MaxwellOperator3D, g::RTRefSpace, f::RTRefSpace, i, τ, j, σ, qd)
  # defines coincidence of points
  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

  # decides on whether to use singularity extraction
  xtol = 0.2

  k = norm(op.gamma)

  hits = 0
  xmin = xtol
  for t in τ.vertices
    for s in σ.vertices
      d = norm(t-s)
      xmin = min(xmin, k*d)
      if d < dtol
        hits +=1
        break
      end
    end
  end

  xmin < xtol && return WiltonSERule(
    qd.tpoints[1,i],
    DoubleQuadRule(
      qd.tpoints[2,i],
      qd.bpoints[2,j],
    ),
  )
  return DoubleQuadRule(
    qd.tpoints[1,i],
    qd.bpoints[1,j],
  )

end
