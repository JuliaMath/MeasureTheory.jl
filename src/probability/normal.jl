
# Normal distribution

import StatsFuns
export Normal

@measure Normal(μ,σ)

basemeasure(::Normal) = (1/sqrt2π) * Lebesgue(ℝ)


basemeasure(d::Normal{(:μ,:σ)}) = ifelse(iszero(d.σ), Dirac(d.μ), (1/sqrt2π) * Lebesgue(ℝ))

basemeasure(d::Normal{(:σ,)}) = ifelse(iszero(d.σ), Dirac(0.0), (1/sqrt2π) * Lebesgue(ℝ))

logdensity(d::Normal{()} , x) = - x^2 / 2 

Base.rand(rng::Random.AbstractRNG, T::Type, μ::Normal{()}) = randn(rng, T)

@μσ_methods Normal()

@half Normal()

HalfNormal(σ) = HalfNormal(σ = σ)

distproxy(d::Normal{(:μ, :σ)}) = Dists.Normal(d.μ, d.σ)
