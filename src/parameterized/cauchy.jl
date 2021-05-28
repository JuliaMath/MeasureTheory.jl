
# Cauchy distribution

import StatsFuns
export Cauchy

@parameterized Cauchy(μ,σ) ≪ (1/π) * Lebesgue(ℝ)

@kwstruct Cauchy()

@μσ_methods Cauchy()
@kwstruct Cauchy(μ,σ)

function logdensity(d::Cauchy{()} , x) 
    return -log(1 + x^2)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Cauchy{()}) = randn(rng, T) / randn(rng, T)

≪(::Cauchy, ::Lebesgue{X}) where X <: Real = true
representative(::Cauchy) = Lebesgue(ℝ)


@half Cauchy()

distproxy(d::Cauchy{()}) = Dists.Cauchy()
