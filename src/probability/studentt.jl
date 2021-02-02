
# StudentT distribution

using SpecialFunctions
using StatsFuns
export StudentT

@measure StudentT(ν) ≃ (1/sqrtπ) * Lebesgue(ℝ)

function logdensity(d::StudentT{(:ν,)}, x) 
    ν = d.ν
    halfνp1 = (ν+1)/2
    return loggamma(halfνp1) - loggamma(ν/2) + ν * log(ν) - halfνp1 * log(x^2 + ν)
end

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

≪(::StudentT, ::Lebesgue{X}) where X <: Real = true
representative(::StudentT) = Lebesgue(ℝ)

@μσ_methods StudentT(ν)
