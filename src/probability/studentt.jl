
# StudentT distribution

using SpecialFunctions
using StatsFuns
export StudentT

import Base: eltype


@measure StudentT(ν) ≃ (1/sqrtπ) * Lebesgue(ℝ)


# Standard StudentT

function logdensity(μ::StudentT{NamedTuple{(:ν,), Tuple{N}}} where {N} , x::X) where {X}  
    ν = μ.par.ν
    halfνp1 = (ν+1)/2
    return loggamma(halfνp1) - loggamma(ν/2) + ν * log(ν) - halfνp1 * log(x^2 + ν)
end

# StudentT() = StudentT{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::StudentT) = Real

Base.rand(μ::StudentT{EmptyNamedTuple}) = rand(Dists.StudentT())

≪(::StudentT, ::Lebesgue{X}) where X <: Real = true
representative(::StudentT) = Lebesgue(ℝ)
