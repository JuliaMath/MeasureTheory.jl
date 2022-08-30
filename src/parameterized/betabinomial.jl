# Beta-Binomial distribution

export BetaBinomial
import Base
using SpecialFunctions

@parameterized BetaBinomial(n, α, β)

basemeasure(d::BetaBinomial) = CountingMeasure()

testvalue(::BetaBinomial) = 0

@kwstruct BetaBinomial(n, α, β)

function Base.rand(rng::AbstractRNG, ::Type, d::BetaBinomial{(:n, :α, :β)})
    rand(rng, Dists.BetaBinomial(d.n, d.α, d.β))
end

function Base.rand(
    rng::AbstractRNG,
    ::Type,
    d::BetaBinomial{(:n, :α, :β),Tuple{I,A}},
) where {I<:Integer,A}
    rand(rng, Dists.BetaBinomial(d.n, d.α, d.β))
end

@inline function insupport(d::BetaBinomial, x)
    isinteger(x) && 0 ≤ x ≤ d.n
end

@inline function logdensity_def(d::BetaBinomial{(:n, :α, :β)}, y)
    (n, α, β) = (d.n, d.α, d.β)
    logbinom = -log1p(n) - logbeta(y + 1, n - y + 1)
    lognum = logbeta(y + α, n - y + β)
    logdenom = logbeta(α, β)
    return logbinom + lognum - logdenom
end

asparams(::Type{<:BetaBinomial}, ::StaticSymbol{:α}) = asℝ₊
asparams(::Type{<:BetaBinomial}, ::StaticSymbol{:β}) = asℝ₊

function proxy(d::BetaBinomial{(:n, :α, :β),Tuple{I,A}}) where {I<:Integer,A}
    Dists.BetaBinomial(d.n, d.α, d.β)
end
