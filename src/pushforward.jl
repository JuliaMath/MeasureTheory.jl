export PushForward

import Bijectors
import Base

struct PushForward{F,M} <: AbstractMeasure
    f::F
    μ::M
end

function logdensity(d::PushForward{F,M}, y) where {F <: Bijectors.Bijector, M}
    res = Bijectors.forward(inv(d.f), y)
    return logdensity(d.μ, res.rv) + res.logabsdetjac
end

for F in [
      Bijectors.ADBijector
    , Bijectors.Composed
    , Bijectors.Exp
    , Bijectors.Identity
    , Bijectors.Inverse
    , Bijectors.InvertibleBatchNorm
    , Bijectors.Log
    , Bijectors.Logit
    , Bijectors.PDBijector
    , Bijectors.Permute
    , Bijectors.PlanarLayer
    , Bijectors.RadialLayer
    , Bijectors.Scale
    , Bijectors.Shift
    , Bijectors.SimplexBijector
    , Bijectors.Stacked
    , Bijectors.TruncatedBijector] 
    @eval (f::$F)(μ::M) where {M <: AbstractMeasure} = PushForward(f,μ)
end



Base.rand(d::PushForward) = d.f(rand(d.μ))
