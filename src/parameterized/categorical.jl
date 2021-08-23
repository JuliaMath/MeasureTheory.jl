# Categorical distribution

# REFERENCES
# https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Categorical
# https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.DiscreteNonParametric

export Categorical

@parameterized Categorical(p) ≪ CountingMeasure(ℤ[0:∞])

ncategories(d::Categorical) = length(d.p)

(d::Categorical ≪ ::CountingMeasure{IntegerRange{a,b}}) where {a,b} = a ≤ 1 && b ≥ ncategories(d)

(::CountingMeasure{IntegerRange{a,b}} ≪ ::Categorical) where {a,b} = a ≥ 1 && b ≤ ncategories(d)

###############################################################################
@kwstruct Categorical(p)

logdensity(d::Categorical{(:p)}, y) = log(d.p[y])

# Very inefficient because of the heavy implementation of Dists.DiscreteNonParametric
distproxy(d::Categorical{(:p)}) = Dists.Categorical(d.p)

Base.rand(rng::AbstractRNG, T::Type, d::Categorical{(:p)}) = rand(rng, distproxy(d))

asparams(::Type{<:Categorical}, ::Val{:p}) = as𝕀
