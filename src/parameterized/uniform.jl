
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


function tbasemeasure_type(::Type{Uniform{(), Tuple{}}}) 
    FactoredBase{typeof(in(𝕀)), Static.StaticFloat64{0.0}, Returns{Static.StaticFloat64{0.0}}, Lebesgue{MeasureBase.RealNumbers}}
end

distproxy(::Uniform{()}) = Dists.Uniform()

logdensity_def(d::Uniform{()}, x) = 0.0

TV.as(::Uniform{()}) = as𝕀

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

###############################################################################
# Uniform
