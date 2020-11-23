
# Exponential distribution

import StatsFuns
export Exponential

import Base: eltype


@measure Exponential(λ) ≪ Lebesgue(ℝ)

# Standard Exponential

function logdensity(d::Exponential{EmptyNamedTuple} , x::X) where {X}  
    return -x
end

# Exponential() = Exponential{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::Exponential{NamedTuple{(),Tuple{}}}) = Real

Base.rand(μ::Exponential{EmptyNamedTuple}) = randn()

≪(::Exponential, ::Lebesgue{X}) where X <: Real = true
representative(::Exponential) = Lebesgue(ℝ)
