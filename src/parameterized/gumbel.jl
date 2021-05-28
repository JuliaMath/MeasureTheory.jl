# Gumbel distribution

import StatsFuns
export Gumbel

@parameterized Gumbel(μ,σ) ≪ Lebesgue(ℝ)

@kwstruct Gumbel()
@kwstruct Gumbel(μ,σ)

@μσ_methods Gumbel()

Gumbel(μ,σ) = Gumbel((μ=μ, σ=σ))

function logdensity(d::Gumbel{()} , x)
    return -exp(-x) - x
end

import Base

function Base.rand(rng::AbstractRNG, d::Gumbel{()})
    u = rand(rng)
    log(-log(u))
end

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true
representative(::Gumbel) = Lebesgue(ℝ)
