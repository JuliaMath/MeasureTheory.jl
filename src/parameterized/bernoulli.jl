# Bernoulli distribution

export Bernoulli
import Base

@parameterized Bernoulli(p) ≃ CountingMeasure(ℤ[0:1])

@inline function logdensity(d::Bernoulli{(:p,)}, y)
    p = d.p
    return y * log(p) + (1 - y) * log(1 - p)
end

function density(d::Bernoulli{(:p,)}, y)
    p = d.p
    return 2 * p * y - p - y + 1
end

@inline function logdensity(d::Bernoulli{(:logitp,)}, y)
    x = d.logitp
    return y * x - log1pexp(x)
end

function density(d::Bernoulli{(:logitp,)}, y)
    exp_x = exp(d.logitp)
    return exp_x^y / (1 + exp_x)
end

sampletype(::Bernoulli) = Bool

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:p,)}) = rand(rng, T) < d.p

Base.rand(rng::AbstractRNG, T::Type, d::Bernoulli{(:logitp,)}) =
    rand(rng, T) < logistic(d.logitp)

≪(::Bernoulli, ::IntegerRange{lo,hi}) where {lo,hi} = lo ≤ 0 && 1 ≤ hi

asparams(::Type{<:Bernoulli}, ::Val{:p}) = as𝕀
asparams(::Type{<:Bernoulli}, ::Val{:logitp}) = asℝ

distproxy(d::Bernoulli{(:p,)}) = Dists.Bernoulli(d.p)
distproxy(d::Bernoulli{(:logitp,)}) = Dists.Bernoulli(logistic(d.logitp))
