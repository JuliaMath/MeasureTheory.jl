
# Normal distribution

import StatsFuns
export Normal

@measure Normal(μ,σ) ≪ (1/sqrt2π) * Lebesgue(ℝ)

logdensity(d::Normal{()} , x) = - x^2 / 2 

Base.rand(rng::Random.AbstractRNG, T::Type, μ::Normal{()}) = randn(rng, T)

@μσ_methods Normal()

@half Normal()
