
# Laplace distribution

import StatsFuns
export Laplace

@parameterized Laplace(μ,σ) ≪ (1/2) * Lebesgue(ℝ)

@kwstruct Laplace()
@kwstruct Laplace(μ,σ)

function logdensity(d::Laplace{()} , x)
    return -abs(x)
end

Base.rand(rng::AbstractRNG, μ::Laplace{()}) = rand(rng, Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true
representative(::Laplace) = Lebesgue(ℝ)

@μσ_methods Laplace()
@half Laplace()
