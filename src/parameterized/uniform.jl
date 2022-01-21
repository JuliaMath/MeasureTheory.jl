
# Uniform distribution

export Uniform

@parameterized Uniform()
@kwstruct Uniform()

###############################################################################
# Standard Uniform


@inline function basemeasure(::Uniform{()})
    inbounds = in(𝕀)
    constℓ = static(0.0)
    varℓ = Returns(static(0.0))
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

distproxy(::Uniform{()}) = Dists.Uniform()

logdensity_def(d::Uniform{()}, x) = 0.0

xform(::Uniform{()}) = as𝕀

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

###############################################################################
# Uniform
