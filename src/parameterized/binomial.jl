# Binomial distribution

export Binomial
import Base
using SpecialFunctions

probit(p) = sqrt2 * erfinv(2p - 1)
Φ(z) = (1 + erf(invsqrt2 * z)) / 2

@parameterized Binomial(n, p) 

basemeasure(d::Binomial) =  Counting(BoundedInts(static(0), d.n))


###############################################################################
@kwstruct Binomial(n, p)

@inline function logdensity_def(d::Binomial{(:n, :p)}, y)
    (n, p) = (d.n, d.p)
    return -log1p(n) - logbeta(n - y + 1, y + 1) + xlogy(y, p) + xlog1py(n - y, -p)
end

function Base.rand(rng::AbstractRNG, ::Type, d::Binomial{(:n, :p)})
    rand(rng, Dists.Binomial(d.n, d.p))
end

###############################################################################
@kwstruct Binomial(n, logitp)

@inline function logdensity_def(d::Binomial{(:n, :logitp)}, y)
    n = d.n
    x = d.logitp
    return -log1p(n) - logbeta(n - y + 1, y + 1) + y * x - n * log1pexp(x)
end

function Base.rand(rng::AbstractRNG, d::Binomial{(:n, :logitp)})
    rand(rng, Dists.Binomial(d.n, logistic(d.logitp)))
end

###############################################################################
@kwstruct Binomial(n, probitp)

@inline function logdensity_def(d::Binomial{(:n, :probitp)}, y)
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

asparams(::Type{<:Binomial}, ::StaticSymbol{:p}) = as𝕀
asparams(::Type{<:Binomial}, ::StaticSymbol{:logitp}) = asℝ
asparams(::Type{<:Binomial}, ::StaticSymbol{:probitp}) = asℝ
