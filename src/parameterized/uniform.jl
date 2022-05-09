
# Uniform distribution

export Uniform

@parameterized Uniform()
@kwstruct Uniform()
@kwstruct Uniform(a, b)

###############################################################################
# Standard Uniform

insupport(::Uniform{()}) = in𝕀
insupport(::Uniform{()}, x) = in𝕀(x)

@inline basemeasure(::Uniform{()}) = Lebesgue(ℝ)

proxy(::Uniform{()}) = Dists.Uniform()

density_def(::Uniform{()}, x) = 1.0

logdensity_def(d::Uniform{()}, x) = 0.0

xform(::Uniform{()}) = as𝕀

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

###############################################################################
# Uniform

@inline insupport(d::Uniform{(:a, :b)}, x) = d.a ≤ x ≤ d.b

Uniform(a, b) = Uniform((a = a, b = b))

proxy(d::Uniform{(:a, :b)}) = affine((μ = d.a, σ = d.b - d.a), Uniform())
@useproxy Uniform{(:a, :b)}
Base.rand(rng::Random.AbstractRNG, ::Type{T}, μ::Uniform) where {T} = rand(rng, T, proxy(μ))
