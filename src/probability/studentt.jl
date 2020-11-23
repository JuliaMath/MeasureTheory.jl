
# StudentT distribution

import SpecialFunctions
export StudentT

import Base: eltype


@measure StudentT(ν) ≃ gamma(0.5) * Lebesgue(ℝ)


# Standard StudentT

function logdensity(d::StudentT{NamedTuple{(:ν,), Tuple{N}}} where {N} , x::X) where {X}  
    νplus1 = ν + 1
    return -0.5 * νplus1 * log(x^2 / νplus1) - log(ν)/2 - loggamma(ν / 2) + loggamma(νplus1 / 2)
end

# StudentT() = StudentT{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::StudentT) = Real

Base.rand(μ::StudentT{EmptyNamedTuple}) = rand(Dists.StudentT())

≪(::StudentT, ::Lebesgue{X}) where X <: Real = true
representative(::StudentT) = Lebesgue(ℝ)
