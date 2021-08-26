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

# The implementation of Dists.DiscreteNonParametric has heavy argument checks
# But I think since the values of Categorical are 1:n the sortperm has no effect
# So it might be OK
distproxy(d::Categorical{(:p)}) = Dists.Categorical(d.p)

Base.rand(rng::AbstractRNG, T::Type, d::Categorical{(:p)}) = rand(rng, distproxy(d))

asparams(::Type{<:Categorical}, ::Val{:p}) = as𝕀

###############################################################################
@kwstruct Categorical(logp)

logdensity(d::Categorical{(:logp)}, y) = d.logp[y]

distproxy(d::Categorical{(:logp)}) = Dists.Categorical(exp.(d.logp))  # inefficient

Base.rand(rng::AbstractRNG, T::Type, d::Categorical{(:logp)}) = rand(rng, distproxy(d))

asparams(::Type{<:Categorical}, ::Val{:logp}) = asℝ
