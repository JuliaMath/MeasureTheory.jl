export PushForward

import Bijectors

const Bij = Bijectors
import Base

struct PushForward{F,M} <: AbstractMeasure
    f::F
    μ::M
end

function logdensity(d::PushForward{F,M}, y) where {F <: Bij.Bijector, M}
    res = Bij.forward(inv(d.f), y)
    return logdensity(d.μ, res.rv) + res.logabsdetjac
end

function logdensity(d::PushForward{F,M}, y) where {F <: TransformVariables.AbstractTransform, M}
    
end

for F in [
      Bij.ADBijector
    , Bij.Composed
    , Bij.Exp
    , Bij.Identity
    , Bij.Inverse
    , Bij.InvertibleBatchNorm
    , Bij.Log
    , Bij.Logit
    , Bij.PDBijector
    , Bij.Permute
    , Bij.PlanarLayer
    , Bij.RadialLayer
    , Bij.Scale
    , Bij.Shift
    , Bij.SimplexBijector
    , Bij.Stacked
    , Bij.TruncatedBijector] 
    @eval (f::$F)(μ::M) where {M <: AbstractMeasure} = PushForward(f,μ)
end



Base.rand(d::PushForward) = d.f(rand(d.μ))
