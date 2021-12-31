# Bernoulli distribution

export Bernoulli
import Base

@parameterized Bernoulli(p)

basemeasure(::Bernoulli) = Counting(Bool)

@inline function logdensity_def(d::Bernoulli{(:p,)}, y)
    p = d.p
    return y * log(p) + (1 - y) * log(1 - p)
end

function density(d::Bernoulli{(:p,)}, y)
    p = d.p
    return 2 * p * y - p - y + 1
end

@inline function logdensity_def(d::Bernoulli{(:logitp,)}, y)
    x = d.logitp
    return y * x - log1pexp(x)
end

function density(d::Bernoulli{(:logitp,)}, y)
    exp_x = exp(d.logitp)
    return exp_x^y / (1 + exp_x)
end

gentype(::Bernoulli) = Bool

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:p,)}) = rand(rng, T) < d.p

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:logitp,)}) =
    rand(rng, T) < logistic(d.logitp)

asparams(::Type{<:Bernoulli}, ::StaticSymbol{:p}) = as𝕀
asparams(::Type{<:Bernoulli}, ::StaticSymbol{:logitp}) = asℝ

distproxy(d::Bernoulli{(:p,)}) = Dists.Bernoulli(d.p)
distproxy(d::Bernoulli{(:logitp,)}) = Dists.Bernoulli(logistic(d.logitp))
