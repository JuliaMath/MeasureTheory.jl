
# Gamma distribution

export Gamma

@parameterized Gamma()

@kwstruct Gamma()

proxy(::Gamma{()}) = Exponential()

@useproxy Gamma{()}

@kwstruct Gamma(k)

@inline function logdensity_def(d::Gamma{(:k,)}, x)
    return xlogy(d.k - 1, x) - x
end

function basemeasure(d::Gamma{(:k,)})
    ℓ = -loggamma(d.k)
    weightedmeasure(ℓ, Lebesgue())
end

@kwstruct Gamma(k, σ)

Gamma(k, σ) = Gamma((k = k, σ = σ))

@useproxy Gamma{(:k, :σ)}
function proxy(d::Gamma{(:k, :σ)})
    affine(NamedTuple{(:σ,)}(d.σ), Gamma((k = d.k,)))
end

@kwstruct Gamma(k, λ)

@useproxy Gamma{(:k, :λ)}
function proxy(d::Gamma{(:k, :λ)})
    affine(NamedTuple{(:λ,)}(d.λ), Gamma((k = d.k,)))
end

Base.rand(rng::AbstractRNG, T::Type, μ::Gamma{()}) = rand(rng, T, Exponential())

Base.rand(rng::AbstractRNG, T::Type, μ::Gamma{(:k,)}) = rand(rng, Dists.Gamma(μ.k))

as(::Gamma) = asℝ₊

insupport(::Gamma, x) = x > 0

@kwstruct Gamma(μ, ϕ)

@inline function logdensity_def(d::Gamma{(:μ, :ϕ)}, x)
    x_μ = x / d.μ
    inv(d.ϕ) * (log(x_μ) - x_μ) - log(x)
end

function basemeasure(d::Gamma{(:μ, :ϕ)})
    ϕ = d.ϕ
    ϕinv = inv(ϕ)
    ℓ = -ϕinv * log(ϕ) - first(logabsgamma(ϕinv))
    weightedmeasure(ℓ, Lebesgue())
end

function basemeasure(d::Gamma{(:μ, :ϕ),Tuple{M,StaticFloat64{ϕ}}}) where {M,ϕ}
    ϕinv = inv(ϕ)
    ℓ = static(-ϕinv * log(ϕ) - first(logabsgamma(ϕinv)))
    weightedmeasure(ℓ, Lebesgue())
end
