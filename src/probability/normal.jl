
# Normal distribution

import StatsFuns
export Normal

@measure Normal(μ,σ) ≃ (1/sqrt2π) * Lebesgue(ℝ)

function logdensity(d::Normal{()} , x::MaybeSym{T}) where {T <: Number}  
    return - x^2 / 2 
end
 
sampletype(::Normal{()}) = Real

Base.rand(rng::Random.AbstractRNG, μ::Normal{()}) = randn(rng)

≪(::Normal, ::Lebesgue{X}) where X <: Real = true
representative(::Normal) = Lebesgue(ℝ)

@μσ_methods Normal()
