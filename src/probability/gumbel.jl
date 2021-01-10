
# Gumbel distribution

import StatsFuns
export Gumbel

import Base: eltype


@measure Gumbel(μ,σ) ≃ Lebesgue(ℝ)


# Standard Gumbel

function logdensity(d::Gumbel{EmptyNamedTuple} , x::X) where {X}  
    return -exp(-x) - x
end

# Gumbel() = Gumbel{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::Gumbel{NamedTuple{(),Tuple{}}}) = Real

Base.rand(μ::Gumbel{EmptyNamedTuple}) = rand(Dists.Gumbel())

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true
representative(::Gumbel) = Lebesgue(ℝ)
