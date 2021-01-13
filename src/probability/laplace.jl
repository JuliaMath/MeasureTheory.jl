
# Laplace distribution

import StatsFuns
export Laplace

import Base: eltype


@measure Laplace(μ,σ) ≃ (1/2) * Lebesgue(ℝ)


# Standard Laplace

function logdensity(d::Laplace{()} , x::X) where {X}  
    return -abs(x)
end

 
sampletype(::Laplace) = Real

Base.rand(rng::AbstractRNG, μ::Laplace{()}) = rand(rng, Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true
representative(::Laplace) = Lebesgue(ℝ)

@μσ_methods Laplace()
