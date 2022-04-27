
# Gamma distribution

export Gamma

@parameterized Gamma(shape) ≃ Lebesgue(ℝ₊)

@inline function logdensity_def(μ::Gamma{(:shape,)}, x)
    α = μ.shape
    xinv = 1 / x

    return xlogy(α + 1, xinv) - xinv - loggamma(α)
end


for N in [(:σ,), (:ω,)]
    @eval begin
        proxy(d::Gamma{(:shape, $N...)}) =
            affine(NamedTuple{$N}(params(d)), Gamma((shape = d.shape,)))
    end
end


Base.rand(rng::AbstractRNG, T::Type, μ::Gamma{(:shape,)}) =
    rand(rng, Dists.Gamma(μ.shape))


TV.as(::Gamma) = asℝ₊

# @μσ_methods Gamma(shape)
