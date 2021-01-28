
# Normal distribution

import StatsFuns
export Normal

@measure Normal(μ,σ) ≪ (1/sqrt2π) * Lebesgue(ℝ)

function logdensity(d::Normal{()} , x)
    return - x^2 / 2 
end

Base.rand(rng::Random.AbstractRNG, T, μ::Normal{()}) = randn(rng, T)

representative(::Normal) = Lebesgue(ℝ)

@μσ_methods Normal()
