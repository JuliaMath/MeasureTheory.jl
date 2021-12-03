
# Cauchy distribution

export Cauchy, HalfCauchy

@parameterized Cauchy(μ, σ)

@kwstruct Cauchy()
@kwstruct Cauchy(μ)
@kwstruct Cauchy(σ)
@kwstruct Cauchy(μ, σ)
@kwstruct Cauchy(ω)
@kwstruct Cauchy(μ, ω)

proxy(d::Cauchy) = affine(params(d), Cauchy())

logdensity_def(d::Cauchy, x) = logdensity_def(proxy(d), x)

basemeasure(d::Cauchy) = (1 / π) * Lebesgue(ℝ)

basemeasure_depth(::Cauchy) = static(2)
basemeasure_depth(::Type{C}) where {C<:Cauchy} = static(2)

# @affinepars Cauchy

@inline function logdensity_def(d::Cauchy{()}, x)
    return -log1p(x^2)
end

function density(d::Cauchy{()}, x)
    return inv(1 + x^2)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Cauchy{()}) = randn(rng, T) / randn(rng, T)

for N in AFFINEPARS
    @eval begin
        rand(rng::AbstractRNG, ::Type{T}, d::Cauchy{$N}) where {T} = rand(rng, T, affine(params(d), Cauchy()))
    end
end

≪(::Cauchy, ::Lebesgue{X}) where {X<:Real} = true

TV.as(::Cauchy) = asℝ

@half Cauchy

HalfCauchy(σ) = HalfCauchy(σ = σ)

distproxy(d::Cauchy{()}) = Dists.Cauchy()
distproxy(d::Cauchy{(:μ,)}) = Dists.Cauchy(d.μ, 1.0)
distproxy(d::Cauchy{(:σ,)}) = Dists.Cauchy(0.0, d.σ)
distproxy(d::Cauchy{(:μ, :σ)}) = Dists.Cauchy(d.μ, d.σ)
distproxy(d::Cauchy{(:ω,)}) = Dists.Cauchy(0.0, inv(d.ω))
distproxy(d::Cauchy{(:μ, :ω)}) = Dists.Cauchy(d.μ, inv(d.ω))
