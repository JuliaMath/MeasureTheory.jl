
# Uniform distribution

export Uniform

@parameterized Uniform()

basemeasure(::Uniform) = Lebesgue(𝕀)

# Standard Uniform

distproxy(::Uniform{()}) = Dists.Uniform()

logdensity(d::Uniform{()}, x) = 0.0

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

representative(::Uniform{()}) = Lebesgue(𝕀)
