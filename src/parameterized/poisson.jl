# Poisson distribution

export Poisson
import Base
using SpecialFunctions: logfactorial

@parameterized Poisson(λ) ≪ CountingMeasure(ℤ[0:∞])

Base.eltype(::Type{P}) where {P<:Poisson} = Int

function logdensity(d::Poisson{(:λ,)}, y)
    λ = d.λ
    return y * log(λ) - λ - logfactorial(y)
end

function logdensity(d::Poisson{(:logλ,)}, y)
    return y * d.logλ + exp(d.logλ) - logfactorial(y)
end

asparams(::Type{<:Poisson}, ::Val{:λ}) = asℝ₊
asparams(::Type{<:Poisson}, ::Val{:logλ}) = asℝ

sampletype(::Poisson) = Int

Base.rand(rng::AbstractRNG, T::Type, d::Poisson{(:λ,)}) = rand(rng, Dists.Poisson(d.λ))
Base.rand(rng::AbstractRNG, T::Type, d::Poisson{(:logλ,)}) = rand(rng, Dists.Poisson(exp(d.logλ)))

≪(::Poisson, ::IntegerRange{lo,hi}) where {lo, hi} = lo ≤ 0 && isinf(hi)
