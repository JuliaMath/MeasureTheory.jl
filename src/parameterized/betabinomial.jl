# Beta-Binomial distribution

export BetaBinomial
import Base
using SpecialFunctions

@parameterized BetaBinomial(n, α, β)

basemeasure(d::BetaBinomial) = CountingBase()

testvalue(::BetaBinomial) = 0

@kwstruct BetaBinomial(n, α, β)

function Base.rand(rng::AbstractRNG, ::Type{T}, d::BetaBinomial{(:n, :α, :β)}) where {T}
    k = rand(rng, T, Beta(d.α, d.β))
    return rand(rng, T, Binomial(d.n, k))
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
