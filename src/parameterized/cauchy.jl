
# Cauchy distribution

export Cauchy, HalfCauchy

@parameterized Cauchy(μ, σ)

@kwstruct Cauchy()
@kwstruct Cauchy(μ)
@kwstruct Cauchy(σ)
@kwstruct Cauchy(μ, σ)
@kwstruct Cauchy(ω)
@kwstruct Cauchy(μ, ω)

logdensity_def(d::Cauchy, x) = logdensity_def(proxy(d), x)

basemeasure(d::Cauchy) = WeightedMeasure(static(-logπ), Lebesgue(ℝ))

# @affinepars Cauchy

@inline function logdensity_def(d::Cauchy{()}, x)
    return -log1p(x^2)
end

function density_def(d::Cauchy{()}, x)
    return inv(1 + x^2)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Cauchy{()}) = randn(rng, T) / randn(rng, T)

for N in AFFINEPARS
    @eval begin
        proxy(d::Cauchy{$N}) = affine(params(d), Cauchy())
        @useproxy Cauchy{$N}
        function rand(rng::AbstractRNG, ::Type{T}, d::Cauchy{$N}) where {T}
            rand(rng, T, affine(params(d), Cauchy()))
        end
    end
end

≪(::Cauchy, ::Lebesgue{X}) where {X<:Real} = true

TV.as(::Cauchy) = asℝ

@half Cauchy

HalfCauchy(σ) = HalfCauchy(σ = σ)

proxy(d::Cauchy{()}) = Dists.Cauchy()

insupport(::Cauchy, x) = true
