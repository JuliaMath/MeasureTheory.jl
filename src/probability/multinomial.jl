# Multinomial distribution

export Multinomial
import Base
using StatsFuns
using SpecialFunctions

probit(p) = sqrt2 * erfinv(2p - 1)
Φ(z) = (1 + erf(invsqrt2 * z))/2

@measure Multinomial(n,p)

basemeasure(μ::Multinomial) = CountingMeasure(ℤ[0:(μ.n)])

(d::Multinomial ≪ ::CountingMeasure{IntegerRange{a,b}}) where {a,b} = a ≤ 0 && b ≥ d.n

(::CountingMeasure{IntegerRange{a,b}} ≪ ::Multinomial) where {a,b} = a ≥ 0 && b ≤ d.n

###############################################################################
# (n, p)
    
function logdensity(d::Multinomial{(:n, :p)}, y)
    (n, p) = (d.n, d.p)
    return -log1p(n) - logbeta(n - y + 1, y + 1) + y * log(p) + (n - y) * log1p(-p)
end

function Base.rand(rng::AbstractRNG, d::Multinomial{(:n,:p)})
    rand(rng, Dists.Multinomial(d.n, d.p))
end

###############################################################################
# (n, logitp)

function logdensity(d::Multinomial{(:n, :logitp)}, y)
    n = d.n
    x = d.logitp
    return  -log1p(n) - logbeta(n - y + 1, y + 1) + y * x - n * log1pexp(x)
end

function Base.rand(rng::AbstractRNG, d::Multinomial{(:n,:logitp)})
    rand(rng, Dists.Multinomial(d.n, logistic(d.logitp)))
end

###############################################################################
# (n, probitp)

function logdensity(d::Multinomial{(:n, :probitp)}, y)
    n = d.n
    z = d.probitp
    return  -log1p(n) - logbeta(n - y + 1, y + 1)  + y * log(Φ(z)) + (n-y) * log(Φ(-z))
end

function Base.rand(rng::AbstractRNG, d::Multinomial{(:n, :probitp)})
    rand(rng, Dists.Multinomial(d.n, Φ(d.probitp)))
end

representative(d::Multinomial) = CountingMeasure(ℤ[0:d.n])

distproxy(d::Multinomial{(:n, :p)}) = Dists.Multinomial(d.n, d.p)
distproxy(d::Multinomial{(:n,:logitp)}) = Dists.Multinomial(d.n, logistic(d.logitp))
distproxy(d::Multinomial{(:n,:probitp)}) = Dists.Multinomial(d.n, Φ(d.probitp))

asparams(::Type{<:Multinomial}, ::Val{:p}) = as𝕀
asparams(::Type{<:Multinomial}, ::Val{:logitp}) = asℝ
asparams(::Type{<:Multinomial}, ::Val{:probitp}) = asℝ
