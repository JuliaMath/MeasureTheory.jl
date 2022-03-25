
# Uniform distribution

export Uniform

@parameterized Uniform()
@kwstruct Uniform()

###############################################################################
# Standard Uniform


insupport(::Uniform{()}) = in𝕀
insupport(::Uniform{()}, x) = in𝕀(x)

@inline function basemeasure(::Uniform{()})
    constℓ = static(0.0)
    varℓ = Returns(static(0.0))
    base = Lebesgue(ℝ)
    FactoredBase(constℓ, varℓ, base)
end
TV.as
distproxy(::Uniform{()}) = Dists.Uniform()

logdensity_def(d::Uniform{()}, x) = 0.0

xform(::Uniform{()}) = as𝕀

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

###############################################################################
# Uniform
