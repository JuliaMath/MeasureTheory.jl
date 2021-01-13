
# Exponential distribution

import StatsFuns
export Exponential

@measure Exponential(λ) ≪ Lebesgue(ℝ)

function logdensity(d::Exponential{EmptyNamedTuple} , x::MaybeSym{T}) where {T <: Number}  
    return -x
end

sampletype(::Exponential{NamedTuple{(),Tuple{}}}) = Real

Base.rand(μ::Exponential{EmptyNamedTuple}) = randn()

≪(::Exponential, ::Lebesgue{X}) where X <: Real = true
representative(::Exponential) = Lebesgue(ℝ)
