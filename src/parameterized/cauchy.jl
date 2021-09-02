
# Cauchy distribution

export Cauchy, HalfCauchy

@parameterized Cauchy(μ,σ) ≃ (1/π) * Lebesgue(ℝ)

@kwstruct Cauchy()

Cauchy(nt::NamedTuple{(:μ,:σ)}) = Affine(nt, Cauchy())
Cauchy(nt::NamedTuple{(:μ,:ω)}) = Affine(nt, Cauchy())
Cauchy(nt::NamedTuple{(:σ,)}) = Affine(nt, Cauchy())
Cauchy(nt::NamedTuple{(:ω,)}) = Affine(nt, Cauchy())
Cauchy(nt::NamedTuple{(:μ,)}) = Affine(nt, Cauchy())

@affinepars Cauchy

function logdensity(d::Cauchy{()} , x) 
    return -log(1 + x^2)
end

function density(d::Cauchy{()} , x) 
    return inv(1 + x^2)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Cauchy{()}) = randn(rng, T) / randn(rng, T)

≪(::Cauchy, ::Lebesgue{X}) where X <: Real = true

TV.as(::Cauchy) = asℝ

@half Cauchy

HalfCauchy(σ) = HalfCauchy(σ=σ)

distproxy(d::Cauchy{()}) = Dists.Cauchy()
