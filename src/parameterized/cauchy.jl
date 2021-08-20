
# Cauchy distribution

export Cauchy

@parameterized Cauchy(μ,σ) ≃ (1/π) * Lebesgue(ℝ)

@kwstruct Cauchy()

@μσ_methods Cauchy()

function logdensity(d::Cauchy{()} , x) 
    return -log(1 + x^2)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Cauchy{()}) = randn(rng, T) / randn(rng, T)

≪(::Cauchy, ::Lebesgue{X}) where X <: Real = true

TV.as(::Cauchy) = asℝ

@half Cauchy()

@kwstruct HalfCauchy()

@σ_methods HalfCauchy()

HalfCauchy(σ) = HalfCauchy(σ=σ)

distproxy(d::Cauchy{()}) = Dists.Cauchy()
