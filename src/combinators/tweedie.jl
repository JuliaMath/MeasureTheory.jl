export Tweedie

abstract type AbstractEDM <: AbstractTransitionKernel end

"""
https://en.wikipedia.org/wiki/Tweedie_distribution

The Tweedie distributions include a number of familiar distributions as well as
some unusual ones, each being specified by the domain of the index parameter. We
have the

    extreme stable distribution, p < 0,
    normal distribution, p = 0,
    Poisson distribution, p = 1,
    compound Poisson-gamma distribution, 1 < p < 2,
    gamma distribution, p = 2,
    positive stable distributions, 2 < p < 3,
    Inverse Gaussian distribution, p = 3,
    positive stable distributions, p > 3, and
    extreme stable distributions, p = ∞.

For 0 < p < 1 no Tweedie model exists. Note that all stable distributions mean
actually generated by stable distributions. 
"""
struct Tweedie{S,B,D,P} <: AbstractEDM
    support_contains::S
    base::B
    dim::D
    p::P
end

struct TweedieMeasure{B,Θ,P,S,C} <: AbstractMeasure
    fam::Tweedie{B,Θ,P}
    θ::Θ
    σ::S
    cumulant::C
end

mean(d::TweedieMeasure) = tweedie_mean(d.fam.p, d.θ)

var(d::TweedieMeasure) = d.σ^2 * mean(d)^d.fam.p

###############################################################################
# Tweedie cumulants

@inline function tweedie_cumulant(p::P, θ) where {P}
    if p == zero(P)
        return 0.5 * θ^2
    elseif p == one(P)
        return exp(θ)
    elseif p == 2
        return -log(-θ)
    else
        α = (p - 2) / (p - 1)
        coeff = (α - 1) / α
        return coeff * (θ / (α - 1))^α
    end
end

@inline function tweedie_cumulant(::StaticFloat64{0.0}, θ)
    return 0.5 * θ^2
end

@inline function tweedie_cumulant(::StaticFloat64{1.0}, θ)
    return exp(θ)
end

@inline function tweedie_cumulant(::StaticFloat64{2.0}, θ)
    return -log(-θ)
end

@generated function tweedie_cumulant(::StaticFloat64{p}, θ) where {p}
    α = (p - 2) / (p - 1)
    coeff = (α - 1) / α

    quote
        $(Expr(:meta, :inline))
        coeff * (θ / (α - 1))^α
    end
end

@inline function (fam::Tweedie)(par)
    base = fam.base(par.σ)
    θ = fam.θ(par)
    η = fam.η(θ)
    t = fam.t
    a = fam.a(θ)
    TweedieMeasure(base, θ, p, σ)
end

###############################################################################
# Tweedie mean function

@inline function tweedie_mean(p::P, θ) where {P}
    if p == zero(P)
        return θ
    elseif p == one(P)
        return exp(θ)
    elseif p == 2
        return inv(log(-θ))
    else
        α_minus_1 = (p - 2) / (p - 1) - 1
        return (θ / α_minus_1)^α_minus_1
    end
end

@inline function tweedie_mean(::StaticFloat64{0.0}, θ)
    return θ
end

@inline function tweedie_mean(::StaticFloat64{1.0}, θ)
    return exp(θ)
end

@inline function tweedie_mean(::StaticFloat64{2.0}, θ)
    return inv(log(-θ))
end

@generated function tweedie_mean(::StaticFloat64{p}, θ) where {p}
    α_minus_1 = (p - 2) / (p - 1) - 1

    quote
        $(Expr(:meta, :inline))
        (θ / α_minus_1)^α_minus_1
    end
end

basemeasure(d::TweedieMeasure) = d.base

function logdensity_def(d::TweedieMeasure, x)
    mydot(x, d.θ) - d.cumulant
end

function MeasureBase.powermeasure(fam::Tweedie, dims::NTuple{N,I}) where {N,I}
    base(σ) = fam.base(σ)^dims
    a = AffineTransform((σ = prod(dims),)) ∘ fam.a
    Tweedie(fam.base^dims, fam.θ, fam.η, t, a)
end

struct TweedieLikelihood{C,Θ,H,T,A} <: AbstractLikelihood
    c::C
    θ::Θ
    η::H
    t::T
    a::A
end

export likelihood

function likelihood(fam::Tweedie, x)
    c = logdensityof(fam.base, x)
    t = fam.t(x)
    TweedieLikelihood(c, fam.θ, fam.η, t, fam.a)
end

@inline function logdensity_def(ℓ::TweedieLikelihood, par)
    θ = ℓ.θ(par)
    mydot(θ, ℓ.t) - ℓ.a(θ) + ℓ.c
end

basemeasure(fam::Tweedie) = fam.base
