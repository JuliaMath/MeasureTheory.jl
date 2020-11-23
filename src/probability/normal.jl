
# Normal distribution

import StatsFuns
export Normal

import Base: eltype


@measure Normal(μ,σ) ≃ (1/sqrt2π) * Lebesgue(ℝ)


function logdensity(d::Normal{P} , x::X) where {P <: NamedTuple{(:μ, :σ)}, X}    
    return - log(d.par.σ)  - (x - d.par.μ)^2 / (2 * d.par.σ^2)
end

# Standard normal

function logdensity(d::Normal{EmptyNamedTuple} , x::X) where {X}  
    return - x^2 / 2 
end

# Normal() = Normal{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::Normal{NamedTuple{(),Tuple{}}}) = Real

Base.rand(μ::Normal{EmptyNamedTuple}) = randn()

≪(::Normal, ::Lebesgue{X}) where X <: Real = true
representative(::Normal) = Lebesgue(ℝ)
