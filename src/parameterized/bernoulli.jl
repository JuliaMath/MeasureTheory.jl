# Bernoulli distribution

export Bernoulli
import Base
using StatsFuns: logistic

@parameterized Bernoulli()

@kwstruct Bernoulli()

@kwstruct Bernoulli(p)

@kwstruct Bernoulli(logitp)

Bernoulli(p) = Bernoulli((p = p,))

basemeasure(::Bernoulli) = CountingBase()

testvalue(::Bernoulli) = true

insupport(::Bernoulli, x) = x == true || x == false

@inline function logdensity_def(d::Bernoulli{(:p,)}, y)
    p = d.p
    y == true ? log(p) : log1p(-p)
end

function density_def(::Bernoulli{()}, y)
    return 0.5
end

@inline function logdensity_def(d::Bernoulli{()}, y)
    return -logtwo
end

function density_def(d::Bernoulli{(:p,)}, y)
    p = d.p
    return 2 * p * y - p - y + 1
end

densityof(d::Bernoulli, y) = density_def(d, y)
unsafe_densityof(d::Bernoulli, y) = density_def(d, y)

@inline function logdensity_def(d::Bernoulli{(:logitp,)}, y)
    x = d.logitp
    return y * x - log1pexp(x)
end

function density_def(d::Bernoulli{(:logitp,)}, y)
    exp_x = exp(d.logitp)
    return exp_x^y / (1 + exp_x)
end

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{()}) = rand(rng, T) < one(T) / 2

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:p,)}) = rand(rng, T) < d.p

function Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:logitp,)})
    rand(rng, T) < logistic(d.logitp)
end

function smf(b::B, x::X) where {B<:Bernoulli,X}
    T = Core.Compiler.return_type(densityof, Tuple{B,X})
    x < zero(x) && return zero(T)
    x ≥ one(x) && return one(T)
    return densityof(b, zero(x))
end

function invsmf(b::Bernoulli, p)
    p0 = densityof(b, 0)
    p ≤ p0 ? 0 : 1
end

proxy(d::Bernoulli{(:p,)}) = Dists.Bernoulli(d.p)
proxy(d::Bernoulli{(:logitp,)}) = Dists.Bernoulli(logistic(d.logitp))
