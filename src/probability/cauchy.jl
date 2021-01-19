
# Cauchy distribution

import StatsFuns
export Cauchy

@measure Cauchy(μ,σ) ≃ (1/π) * Lebesgue(ℝ)

function logdensity(d::Cauchy{()} , x) 
    return -log(1 + x^2)
end

Base.rand(μ::Cauchy{()}) = rand(Dists.Cauchy())

≪(::Cauchy, ::Lebesgue{X}) where X <: Real = true
representative(::Cauchy) = Lebesgue(ℝ)

@μσ_methods Cauchy()
