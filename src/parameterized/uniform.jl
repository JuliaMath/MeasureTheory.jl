
# Uniform distribution

export Uniform

@parameterized Uniform()
@kwstruct Uniform()


###############################################################################
# Standard Uniform

function basemeasure(::Uniform{()})
    inbounds(x) = 0 < x < 1
    constℓ = 0.0
    varℓ = 0.0
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

distproxy(::Uniform{()}) = Dists.Uniform()

logdensity(d::Uniform{()}, x) = 0.0

TV.as(::Uniform{()}) = as𝕀

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

###############################################################################
# Uniform
