
# Cauchy distribution

import StatsFuns
export Cauchy

@measure Cauchy(μ,σ) 

basemeasure(::Cauchy) = (1/π) * Lebesgue(ℝ)

function logdensity(d::Cauchy{()} , x) 
    return -log(1 + x^2)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Cauchy{()}) = randn(rng, T) / randn(rng, T)

≪(::Cauchy, ::Lebesgue{X}) where X <: Real = true
representative(::Cauchy) = Lebesgue(ℝ)

@μσ_methods Cauchy()
@half Cauchy()

distproxy(d::Cauchy{()}) = Dists.Cauchy()
