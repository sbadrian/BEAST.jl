
struct HH2DDoubleLayerNear{K}
  gamma::K
end

struct HH2DSingleLayerNear{K}
  gamma::K
end

struct HH2DHyperSingularNear{K}
  gamma::K
end

struct HH2DDoubleLayerTransposedNear{K}
  gamma::K
end

function HH2DDoubleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DDoubleLayerNear(-imag(wavenumber))
  else
    HH2DDoubleLayerNear(wavenumber*im)
  end
end

function HH2DSingleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DSingleLayerNear(-imag(wavenumber))
  else
    HH2DSingleLayerNear(wavenumber*im)
  end
end

function HH2DHyperSingularNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DHyperSingularNear(-imag(wavenumber))
  else
    HH2DHyperSingularNear(wavenumber*im)
  end
end

function HH2DDoubleLayerTransposedNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH2DDoubleLayerTransposedNear(-imag(wavenumber))
  else
    HH2DDoubleLayerTranposedNear(wavenumber*im)
  end
end

HH2DNear = Union{HH2DSingleLayerNear, HH2DDoubleLayerNear, HH2DDoubleLayerTransposedNear, HH2DHyperSingularNear}
defaultquadstrat(op::HH2DNear, basis) = SingleNumQStrat(20)
quaddata(op::HH2DNear,rs,els,qs::SingleNumQStrat) = quadpoints(rs,els,(qs.quad_rule,))
quadrule(op::HH2DNear,refspace,p,y,q,el,qdata,qs::SingleNumQStrat) = qdata[1,q]

function kernelvals(op::HH2DNear,y,p)

    if iszero(real(op.gamma))
        k = imag(op.gamma)
    else
        k = -im*op.gamma
    end
    x = cartesian(p)
    r = y - x
    R = norm(r)

    kr = k * R
    hankels = hankelh2.([0 1], kr)
    #green = - op.alpha * im / 4 * hankels[1]
    #gradgreen = op.alpha * k * im / 4 * hankels[2] * r / R
    green = - im / 4 * hankels[1]
    gradgreen =  k * im / 4 * hankels[2] * r / R


    #txty = dot(normal(tgeo), normal(bgeo))

    nx = normal(p)

    (;op.gamma, r, R, green, gradgreen, nx)
end

function integrand(op::HH2DDoubleLayerNear,krn,y,f,p)

  ∇G = -krn.gradgreen
  nx = krn.nx

  fx = f.value

  ∂G∂n = nx ⋅ ∇G

  return ∂G∂n * fx 
end

function integrand(op::HH2DSingleLayerNear, krn, y, f, p)
  G = krn.green
  fx = f.value

  return G * fx
end

function integrand(op::HH2DDoubleLayerTransposedNear, krn, y, f, p)
  ∇G = krn.gradgreen
  fx = f.value

  return ∇G * fx
end

function integrand(op::HH2DHyperSingularNear, krn, y, f, p)
  G = krn.green
  nx = krn.nx
  γ = krn.γ
  invR = 1/krn.R
  fx = f.value
  r = krn.r

  # returns ∇ᵣn̂'∇ᵣ'G - strong singularity, probably makes problems for small distances R
  return (nx * G * invR * (γ + invR) - r * dot(nx,r) * G * invR^2 * (γ^2 + 3 * γ * invR + 3 * invR^2)) * fx
end
