
# Snedecor's F distribution

export SnedecorF

@parameterized SnedecorF(ν1, ν2)

@kwstruct SnedecorF(ν1, ν2)

@inline function logdensity_def(d::SnedecorF{(:ν1, :ν2)}, y)
    ν1, ν2 = d.ν1, d.ν2
    ν1ν2 = ν1 / ν2
    val = (xlogy(ν1, ν1ν2) + xlogy(ν1 - 2, y) - xlogy(ν1 + ν2, 1 + ν1ν2 * y)) / 2
    return val
end

@inline function basemeasure(d::SnedecorF{(:ν1, :ν2)})
    ℓ = -logbeta(d.ν1 / 2, d.ν2 / 2)
    weightedmeasure(ℓ, Lebesgue())
end

xform(::SnedecorF) = asℝ₊

Base.rand(rng::AbstractRNG, T::Type, d::SnedecorF) = rand(rng, proxy(d))

proxy(d::SnedecorF{(:ν1, :ν2)}) = Dists.FDist(d.ν1, d.ν2)

asparams(::Type{<:SnedecorF}, ::StaticSymbol{:ν1}) = asℝ₊
asparams(::Type{<:SnedecorF}, ::StaticSymbol{:ν2}) = asℝ₊

insupport(::SnedecorF, x) = x > 0

# cdf(d::SnedecorF, x) = StatsFuns.fdistcdf(d.ν1, d.ν2, x)
# ccdf(d::SnedecorF, x) = StatsFuns.fdistccdf(d.ν1, d.ν2, x)