
# Normal distribution

import StatsFuns
export Normal

import Base: eltype


@measure Normal(μ,σ) ≃ (1/sqrt2π) * Lebesgue(Real)


function logdensity(d::Normal{(:μ, :σ)} , x::X) where {X}    
    return - log(d.σ)  - (x - d.μ)^2 / (2 * d.σ^2)
end

function Base.rand(rng::Random.AbstractRNG, d::Normal{(:μ, :σ)})   
    return randn(rng) * d.σ + d.μ
end

# Standard normal

function logdensity(d::Normal{()} , x::X) where {X}  
    return - x^2 / 2 
end

# Normal() = Normal{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::Normal{()}) = Real

Base.rand(rng::Random.AbstractRNG, μ::Normal{()}) = randn(rng)

≪(::Normal, ::Lebesgue{X}) where X <: Real = true
representative(::Normal) = Lebesgue(Real)
