export Tweedie

struct Tweedie{B,Θ,P} <: AbstractKleisli
    base::B
    θ::Θ
    p::P
end

struct TweedieMeasure{B,Θ,P,S,C} <: AbstractMeasure
    base::B
    θ::Θ
    p::P
    σ::S
    cumulant::C
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
