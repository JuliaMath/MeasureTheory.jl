# Binomial distribution

export Binomial
import Base
using SpecialFunctions

probit(p) = sqrt2 * erfinv(2p - 1)
Φ(z) = (1 + erf(invsqrt2 * z)) / 2

@parameterized Binomial(n, p) ≪ CountingMeasure(ℤ[0:∞])

(d::Binomial ≪ ::CountingMeasure{IntegerRange{a,b}}) where {a,b} = a ≤ 0 && b ≥ d.n

(::CountingMeasure{IntegerRange{a,b}} ≪ ::Binomial) where {a,b} = a ≥ 0 && b ≤ d.n

###############################################################################
@kwstruct Binomial(n, p)

@inline function logdensity(d::Binomial{(:n, :p)}, y)
    (n, p) = (d.n, d.p)
    return -log1p(n) - logbeta(n - y + 1, y + 1) + xlogy(y, p) + xlog1py(n - y, -p)
end

function Base.rand(rng::AbstractRNG, T::Type, d::Binomial{(:n, :p)})
    rand(rng, Dists.Binomial(d.n, d.p))
end

###############################################################################
@kwstruct Binomial(n, logitp)

@inline function logdensity(d::Binomial{(:n, :logitp)}, y)
    n = d.n
    x = d.logitp
    return -log1p(n) - logbeta(n - y + 1, y + 1) + y * x - n * log1pexp(x)
end

function Base.rand(rng::AbstractRNG, d::Binomial{(:n, :logitp)})
    rand(rng, Dists.Binomial(d.n, logistic(d.logitp)))
end

###############################################################################
@kwstruct Binomial(n, probitp)

@inline function logdensity(d::Binomial{(:n, :probitp)}, y)
    n = d.n
    z = d.probitp
    return -log1p(n) - logbeta(n - y + 1, y + 1) + xlogy(y, Φ(z)) + xlogy(n - y, Φ(-z))
end

function Base.rand(rng::AbstractRNG, d::Binomial{(:n, :probitp)})
    rand(rng, Dists.Binomial(d.n, Φ(d.probitp)))
end

distproxy(d::Binomial{(:n, :p)}) = Dists.Binomial(d.n, d.p)
distproxy(d::Binomial{(:n, :logitp)}) = Dists.Binomial(d.n, logistic(d.logitp))
distproxy(d::Binomial{(:n, :probitp)}) = Dists.Binomial(d.n, Φ(d.probitp))

asparams(::Type{<:Binomial}, ::Val{:p}) = as𝕀
asparams(::Type{<:Binomial}, ::Val{:logitp}) = asℝ
asparams(::Type{<:Binomial}, ::Val{:probitp}) = asℝ
