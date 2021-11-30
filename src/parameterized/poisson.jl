# Poisson distribution

export Poisson
import Base
using SpecialFunctions: logfactorial

@parameterized Poisson(λ) ≪ CountingMeasure(ℤ[0:∞])

Base.eltype(::Type{P}) where {P<:Poisson} = Int

@inline function logdensity_def(d::Poisson{(:λ,)}, y)
    λ = d.λ
    return y * log(λ) - λ - logfactorial(y)
end

@inline function logdensity_def(d::Poisson{(:logλ,)}, y)
    return y * d.logλ - exp(d.logλ) - logfactorial(y)
end

asparams(::Type{<:Poisson}, ::StaticSymbol{:λ}) = asℝ₊
asparams(::Type{<:Poisson}, ::StaticSymbol{:logλ}) = asℝ

gentype(::Poisson) = Int

Base.rand(rng::AbstractRNG, T::Type, d::Poisson{(:λ,)}) = rand(rng, Dists.Poisson(d.λ))
Base.rand(rng::AbstractRNG, T::Type, d::Poisson{(:logλ,)}) =
    rand(rng, Dists.Poisson(exp(d.logλ)))
