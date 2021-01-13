# Gumbel distribution

import StatsFuns
export Gumbel

@measure Gumbel(μ,σ) ≃ Lebesgue(ℝ)

function logdensity(d::Gumbel{EmptyNamedTuple} , x::MaybeSym{T}) where {T <: Number}  
    return -exp(-x) - x
end
 
sampletype(::Gumbel{NamedTuple{(),Tuple{}}}) = Real

Base.rand(μ::Gumbel{EmptyNamedTuple}) = rand(Dists.Gumbel())

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true
representative(::Gumbel) = Lebesgue(ℝ)
