
# Uniform distribution

export Uniform

@measure Uniform() ≃ Lebesgue(𝕀)

# Standard Uniform

logdensity(d::Uniform{()}, x) = 0.0

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

representative(::Uniform{()}) = Lebesgue(𝕀)
