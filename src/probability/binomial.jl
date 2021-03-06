# Binomial distribution

export Binomial
import Base
using StatsFuns
using SpecialFunctions

probit(p) = sqrt2 * erfinv(2p - 1)
Φ(z) = (1 + erf(invsqrt2 * z))/2

struct Binomial{N, T} <: ParameterizedMeasure{N, T}
    par::NamedTuple{N, T}
end

Binomial(nt::NamedTuple{N,T}) where {N, T} = Binomial{N,T}(nt)

basemeasure(μ::Binomial) = CountingMeasure(ℤ[0:(μ.n)])

Binomial(n,p) = Binomial(; n, p)

(d::Binomial ≪ ::CountingMeasure{IntegerRange{a,b}}) where {a,b} = a ≤ 0 && b ≥ d.n

(::CountingMeasure{IntegerRange{a,b}} ≪ ::Binomial) where {a,b} = a ≥ 0 && b ≤ d.n

###############################################################################
# (n, p)
    
function logdensity(d::Binomial{(:n, :p)}, y)
    (n, p) = (d.n, d.p)
    return -log1p(n) - logbeta(n - y + 1, y + 1) + y * log(p) + (n - y) * log1p(-p)
end

function Base.rand(rng::AbstractRNG, T::Type, d::Binomial{(:n,:p)})
    rand(rng, Dists.Binomial(d.n, d.p))
end

###############################################################################
# (n, logit_p)

function logdensity(d::Binomial{(:n, :logit_p)}, y)
    n = d.n
    x = d.logit_p
    return  -log1p(n) - logbeta(n - y + 1, y + 1) + y * x - n * log1pexp(x)
end

function Base.rand(rng::AbstractRNG, T::Type, d::Binomial{(:n,:logit_p)})
    rand(rng, Dists.Binomial(d.n, logistic(d.logit_p)))
end

###############################################################################
# (n, probit_p)

function logdensity(d::Binomial{(:n, :probit_p)}, y)
    n = d.n
    z = d.probit_p
    return  -log1p(n) - logbeta(n - y + 1, y + 1)  + y * log(Φ(z)) + (n-y) * log(Φ(-z))
end

function Base.rand(rng::AbstractRNG, T::Type, d::Binomial{(:n, :probit_p)})
    rand(rng, Dists.Binomial(d.n, Φ(d.probit_p)))
end

representative(d::Binomial) = CountingMeasure(ℤ[0:d.n])

distproxy(d::Binomial{(:n, :p)}) = Dists.Binomial(d.n, d.p)
distproxy(d::Binomial{(:n,:logit_p)}) = Dists.Binomial(d.n, logistic(d.logit_p))
distproxy(d::Binomial{(:n,:probit_p)}) = Dists.Binomial(d.n, Φ(d.probit_p))
