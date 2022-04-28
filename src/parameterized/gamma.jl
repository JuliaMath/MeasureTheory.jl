
# Gamma distribution

export Gamma

@parameterized Gamma()

@kwstruct Gamma()

@inline function logdensity_def(d::Gamma{()}, x)
    return -x
end

@kwstruct Gamma(k)

@inline function logdensity_def(d::Gamma{(:k,)}, x)
    return xlogy(d.k - 1, x) - x
end

function basemeasure(d::Gamma{(:k,)})
    constℓ = 0.0
    varℓ() = -loggamma(d.k)
    base = LebesgueMeasure()
    FactoredBase(constℓ, varℓ, base)
end

@kwstruct Gamma(k, θ)

Gamma(k,θ) = Gamma((k=k, θ=θ))

@useproxy Gamma{(:k, :θ)}
function proxy(d::Gamma{(:k, :θ)})
    affine(NamedTuple{(:σ,)}(d.θ), Gamma((k = d.k,)))
end

Base.rand(rng::AbstractRNG, T::Type, μ::Gamma{()}) =
    rand(rng, T, Exponential())

Base.rand(rng::AbstractRNG, T::Type, μ::Gamma{(:k,)}) =
    rand(rng, Dists.Gamma(μ.k))


TV.as(::Gamma) = asℝ₊

insupport(::Gamma, x) = x > 0

@kwstruct Gamma(μ, ϕ)

@inline function logdensity_def(d::Gamma{(:μ,:ϕ)}, x)
    x_μ = x / d.μ
    inv(d.ϕ) * (log(x_μ) - x_μ) - log(x)
end

function basemeasure(d::Gamma{(:μ,:ϕ)})
    constℓ = 0.0
    function varℓ()
        ϕ = dynamic(d.ϕ)
        ϕinv = inv(ϕ)
        -ϕinv * log(ϕ) - first(logabsgamma(ϕinv))
    end
    base = LebesgueMeasure()
    FactoredBase(constℓ, varℓ, base)
end