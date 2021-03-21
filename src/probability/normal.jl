
# Normal distribution

import StatsFuns
export Normal

@measure Normal(μ,σ)

basemeasure(::Normal) = (1/sqrt2π) * Lebesgue(ℝ)

logdensity(d::Normal{()} , x) = - x^2 / 2 

Base.rand(rng::Random.AbstractRNG, T::Type, μ::Normal{()}) = randn(rng, T)

@μσ_methods Normal()

@half Normal()

HalfNormal(σ) = HalfNormal(σ = σ)

distproxy(d::Normal{(:μ, :σ)}) = Dists.Normal(d.μ, d.σ)

asparams(::Type{<:Normal}, ::Val{:μ}) = asℝ
asparams(::Type{<:Normal}, ::Val{:σ}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:logσ}) = asℝ
asparams(::Type{<:Normal}, ::Val{:σ²}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:τ}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:logτ}) = asℝ
