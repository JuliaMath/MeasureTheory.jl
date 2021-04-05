
# StudentT distribution

using SpecialFunctions
using StatsFuns
export StudentT

@measure StudentT(ν)

basemeasure(::StudentT) = (1/sqrtπ) * Lebesgue(ℝ)

function logdensity(d::StudentT{(:ν,)}, x) 
    ν = d.ν
    halfνp1 = (ν+1)/2
    return loggamma(halfνp1) - loggamma(ν/2) + ν * log(ν) - halfνp1 * log(x^2 + ν)
end

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, T, Dists.TDist(μ.ν))

≪(::StudentT, ::Lebesgue{X}) where X <: Real = true
representative(::StudentT) = Lebesgue(ℝ)

distproxy(d::StudentT{(:ν, :μ, :σ)}) = Dists.LocationScale(d.μ, d.σ, Dists.TDist(d.ν))

@μσ_methods StudentT(ν)
@half StudentT()
