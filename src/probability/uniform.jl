
# Uniform distribution

export Uniform

@measure Uniform()

basemeasure(::Uniform) = Lebesgue(𝕀)

# Standard Uniform

logdensity(d::Uniform{()}, x) = 0.0

Base.rand(rng::AbstractRNG, μ::Uniform{()}) = rand(rng)

representative(::Uniform{()}) = Lebesgue(𝕀)
