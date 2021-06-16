# Gumbel distribution

import StatsFuns
export Gumbel

@parameterized Gumbel(μ,σ) ≃ Lebesgue(ℝ)

@kwstruct Gumbel()

@μσ_methods Gumbel()

function logdensity(d::Gumbel{()} , x)
    return -exp(-x) - x
end

import Base

function Base.rand(rng::AbstractRNG, d::Gumbel{()})
    u = rand(rng)
    -log(-log(u))
end

TV.as(::Gumbel) = asℝ

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true
