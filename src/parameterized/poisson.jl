# Poisson distribution

export Poisson
import Base
using SpecialFunctions: loggamma

@parameterized Poisson()

@kwstruct Poisson()
@kwstruct Poisson(λ)

Poisson(λ) = Poisson((λ = λ,))

basemeasure(::Poisson) = Counting(BoundedInts(static(0), static(Inf)))

@inline function logdensity_def(d::Poisson{()}, y::T) where {T}
    return y - one(T) - loggamma(one(T) + y)
end

@inline function logdensity_def(d::Poisson{(:λ,)}, y)
    λ = d.λ
    return y * log(λ) - λ - loggamma(1 + y)
end

@kwstruct Poisson(logλ)

@inline function logdensity_def(d::Poisson{(:logλ,)}, y)
    return y * d.logλ - exp(d.logλ) - loggamma(1 + y)
end

gentype(::Poisson) = Int

Base.rand(rng::AbstractRNG, T::Type, d::Poisson{(:λ,)}) = rand(rng, Dists.Poisson(d.λ))
function Base.rand(rng::AbstractRNG, T::Type, d::Poisson{(:logλ,)})
    rand(rng, Dists.Poisson(exp(d.logλ)))
end

mean(d::Poisson{(:λ,)}) = d.λ
std(d::Poisson{(:λ,)}) = sqrt(d.λ)
var(d::Poisson{(:λ,)}) = d.λ

@inline function insupport(::Poisson, x)
    isinteger(x) && x ≥ 0
end
