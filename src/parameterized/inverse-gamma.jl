
# InverseGamma distribution

export InverseGamma

@parameterized InverseGamma(shape)

@kwstruct InverseGamma(shape)

@kwstruct InverseGamma(shape, σ)

insupport(::InverseGamma, x::Real) = x > 0

basemeasure(::InverseGamma) = LebesgueBase()

@inline function logdensity_def(μ::InverseGamma{(:shape,)}, x)
    α = μ.shape
    xinv = 1 / x

    α - loggamma(α) + xlogy(α + 1, xinv) - xinv
end

function Base.rand(rng::AbstractRNG, T::Type, μ::InverseGamma{(:shape,)})
    rand(rng, Dists.InverseGamma(μ.shape, one(μ.shape)))
end

InverseGamma(shape, σ) = InverseGamma((shape = shape, σ = σ))

proxy(d::InverseGamma{(:shape, :σ)}) = affine((σ = d.σ,), InverseGamma(d.shape))

function Base.rand(rng::AbstractRNG, T::Type, d::InverseGamma{(:shape, :σ)})
    rand(rng, T, proxy(d))
end

@useproxy InverseGamma{(:shape, :σ)}
