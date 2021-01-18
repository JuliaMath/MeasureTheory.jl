# Gumbel distribution

import StatsFuns
export Gumbel

@measure Gumbel(μ,σ) ≃ Lebesgue(ℝ)

function logdensity(d::Gumbel{()} , x)
    return -exp(-x) - x
end
 
sampletype(::Gumbel{()}) = Real

import Base

function Base.rand(rng::AbstractRNG, d::Gumbel{()})
    u = rand(rng)
    log(-log(u))
end

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true
representative(::Gumbel) = Lebesgue(ℝ)

@μσ_methods Gumbel()
