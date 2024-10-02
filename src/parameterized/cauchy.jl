
# Cauchy distribution

export Cauchy, HalfCauchy

@parameterized Cauchy(μ, σ)

@kwstruct Cauchy()
@kwstruct Cauchy(μ)
@kwstruct Cauchy(σ)
@kwstruct Cauchy(μ, σ)
@kwstruct Cauchy(λ)
@kwstruct Cauchy(μ, λ)

logdensity_def(d::Cauchy, x) = logdensity_def(proxy(d), x)

const logπ = Float64(log(big(π)))

basemeasure(d::Cauchy) = WeightedMeasure(static(-logπ), LebesgueBase())

# @affinepars Cauchy

@inline function logdensity_def(d::Cauchy{()}, x)
    return -log1p(x^2)
end

function density_def(d::Cauchy{()}, x)
    return inv(1 + x^2)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Cauchy{()}) = randn(rng, T) / randn(rng, T)

Base.rand(::ConstantRNG, ::Type{T}, ::Cauchy{()}) where {T} = zero(T)

for N in AFFINEPARS
    @eval begin
        proxy(d::Cauchy{$N}) = affine(params(d), Cauchy())
        @useproxy Cauchy{$N}
        function rand(rng::AbstractRNG, ::Type{T}, d::Cauchy{$N}) where {T}
            rand(rng, T, affine(params(d), Cauchy()))
        end
    end
end


@half Cauchy

HalfCauchy(σ) = HalfCauchy(σ = σ)

Base.rand(::ConstantRNG, ::Type{T}, μ::Half{<:Cauchy}) where {T} = one(T)

insupport(::Cauchy, x) = true

function smf(::Cauchy{()}, x)
    invπ * atan(x) + 0.5
end

function invsmf(::Cauchy{()}, p)
    tan(π * (p - 0.5))
end

proxy(d::Cauchy{()}) = Dists.Cauchy()
