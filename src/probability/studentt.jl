
# StudentT distribution

using SpecialFunctions
using StatsFuns
export StudentT

@measure StudentT(ν) ≃ (1/sqrtπ) * Lebesgue(ℝ)

function logdensity(μ::StudentT{(:ν,)}, x::MaybeSym{T}) where {T <: Number}  
    ν = μ.ν
    halfνp1 = (ν+1)/2
    return loggamma(halfνp1) - loggamma(ν/2) + ν * log(ν) - halfνp1 * log(x^2 + ν)
end

sampletype(::StudentT) = Real

Base.rand(rng::AbstractRNG, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

≪(::StudentT, ::Lebesgue{X}) where X <: Real = true
representative(::StudentT) = Lebesgue(ℝ)

@μσ_methods StudentT(ν)
