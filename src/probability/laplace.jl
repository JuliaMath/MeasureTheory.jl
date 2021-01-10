
# Laplace distribution

import StatsFuns
export Laplace

import Base: eltype


@measure Laplace(μ,σ) ≃ (1/2) * Lebesgue(ℝ)


# Standard Laplace

function logdensity(d::Laplace{EmptyNamedTuple} , x::X) where {X}  
    return -abs(x)
end

# Laplace() = Laplace{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::Laplace{NamedTuple{(),Tuple{}}}) = Real

Base.rand(μ::Laplace{EmptyNamedTuple}) = rand(Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true
representative(::Laplace) = Lebesgue(ℝ)
