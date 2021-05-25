# Poisson distribution

export Poisson
import Base
using SpecialFunctions: logfactorial

@measure Poisson(λ)

basemeasure(::Poisson) = CountingMeasure(ℤ[0:∞])

function logdensity(d::Poisson{(:λ,)}, y)
    λ = d.λ
    return y * log(λ) - λ - logfactorial(y)
end

# function logdensity(d::Poisson{(:log_λ,)}, y)
#     log_λ = d.log_λ
#     return y * log(log_λ) + (1 - y) * log(1 - log_λ)
# end

# function density(d::Poisson{(:log_λ,)}, y)
#     log_λ = d.log_λ
#     return 2*log_λ*y - log_λ - y + 1
# end

sampletype(::Poisson) = Int

Base.rand(rng::AbstractRNG, T::Type, d::Poisson{(:λ,)}) = rand(rng, T, Dists.Poisson(d.λ))

≪(::Poisson, ::IntegerRange{lo,hi}) where {lo, hi} = lo ≤ 0 && isinf(hi)

representative(::Poisson) = CountingMeasure(ℤ[0:∞])
