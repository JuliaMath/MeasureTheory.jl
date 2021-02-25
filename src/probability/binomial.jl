# Binomial distribution

export Binomial
import Base
using StatsFuns
using SpecialFunctions

probit(p) = sqrt2 * erfinv(2p - 1)
Φ(z) = (1 + erf(invsqrt2 * z))/2

struct Binomial{n, N, T} <: ParameterizedMeasure{N, T}
    par::NamedTuple{N, T}
end

Binomial{n}(nt::NamedTuple{N,T}) where {n, N, T} = Binomial{n,N,T}(nt)

function MeasureTheory.basemeasure(μ::Binomial{n}) where {n}
    return WeightedMeasure(-log1p(n) + loggamma(n + 2), CountingMeasure(ℤ[0:n]))
end

(::Binomial{n} ≪ ::CountingMeasure{IntegerRange{a,b}}) where {n,a,b} = a ≤ 0 && b ≥ n

(::CountingMeasure{IntegerRange{a,b}} ≪ ::Binomial{n}) where {n,a,b} = a ≥ 0 && b ≤ n

Binomial(n,p) = Binomial{n}(; n, p)


Binomial(nt::NamedTuple) = Binomial{nt.n}(nt)
    
function logdensity(d::Binomial{n, (:n, :p)}, y) where {n}
    p = d.p
    return  - loggamma(n - y + 1) - loggamma(y + 1)  + y * log(p) + (n - y) * log1p(-p)
end

function logdensity(d::Binomial{n, (:n, :logit_p)}, y) where {n}
    x = d.logit_p
    return  - loggamma(n - y + 1) - loggamma(y + 1)  + y * x - n * log1pexp(x)
end

function logdensity(d::Binomial{n, (:n, :probit_p)}, y) where {n}
    z = d.probit_p
    return  - loggamma(n - y + 1) - loggamma(y + 1)  + y * log(Φ(z)) + (n-y) * log(Φ(-z))
end

sampletype(::Binomial) = Int

function Base.rand(rng::AbstractRNG, T::Type, d::Binomial{n, (:n,:p)}) where {n}
    rand(rng, Dists.Binomial(d.n, d.p))
end

function Base.rand(rng::AbstractRNG, T::Type, d::Binomial{n, (:n,:logit_p)}) where {n}
    rand(rng, Dists.Binomial(d.n, logistic(d.logit_p)))
end

function Base.rand(rng::AbstractRNG, T::Type, d::Binomial{n, (:n,:probit_p)}) where {n}
    rand(rng, Dists.Binomial(d.n, Φ(d.probit_p)))
end

representative(::Binomial{n}) where {n} = CountingMeasure(ℤ[0:n])


distproxy(d::Binomial{n, (:n, :p)}) where {n} = Dists.Binomial(d.n, d.p)
