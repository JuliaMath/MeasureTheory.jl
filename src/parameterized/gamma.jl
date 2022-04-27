
# Gamma distribution

export Gamma

@parameterized Gamma(k)

@kwstruct Gamma(k)
@kwstruct Gamma(k, θ)

Gamma(k,θ) = Gamma((k=k, θ=θ))

@useproxy Gamma{(:k, :θ)}

@inline function logdensity_def(d::Gamma{(:k,)}, x)
    return xlogy(d.k - 1, x) - x
end

function basemeasure(d::Gamma{(:k,)})
    constℓ = 0.0
    varℓ() = -loggamma(d.k)
    base = LebesgueMeasure()
    FactoredBase(constℓ, varℓ, base)
end

function proxy(d::Gamma{(:k, :θ)})
    affine(NamedTuple{(:σ,)}(d.θ), Gamma((k = d.k,)))
end


Base.rand(rng::AbstractRNG, T::Type, μ::Gamma{(:k,)}) =
    rand(rng, Dists.Gamma(μ.k))


TV.as(::Gamma) = asℝ₊

insupport(::Gamma, x) = x > 0
