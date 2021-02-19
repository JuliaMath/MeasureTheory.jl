# Bernoulli distribution

export Bernoulli
import Base

@measure Bernoulli(p) ≃ CountingMeasure(ℤ[0:1])

function logdensity(d::Bernoulli{(:p,)}, y)
    p = d.p
    return y * log(p) + (1 - y) * log(1 - p)
end

function density(d::Bernoulli{(:p,)}, y)
    p = d.p
    return 2*p*y - p - y + 1
end

function logdensity(d::Bernoulli{(:logit_p,)}, y)
    x = d.logit_p
    return y * x - log1pexp(x)
end

function density(d::Bernoulli{(:logit_p,)}, y)
    exp_x = exp(d.logit_p)
    return exp_x ^ y / (1 + exp_x)
end

sampletype(::Bernoulli) = Bool

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:p,)}) = T(rand(rng) < d.p)

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:logit_p,)}) = T(rand(rng) < logistic(d.logit_p))

≪(::Bernoulli, ::IntegerRange{lo,hi}) where {lo, hi} = lo ≤ 0 && 1 ≤ hi

representative(::Bernoulli) = CountingMeasure(ℤ[0:1])

distproxy(d::Bernoulli{(:p,)}) = Dists.Bernoulli(d.p)
distproxy(d::Bernoulli{(:logit_p,)}) = Dists.Bernoulli(logistic(d.logit_p))
