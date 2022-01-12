export ExponentialFamily

struct ExponentialFamily{B,Θ,H,T,A} <: AbstractKleisli
    base::B
    θ::Θ
    η::H
    t::T
    a::A
end

struct ExpFamMeasure{B,Θ,H,T,A,P} <: AbstractMeasure
    base::B
    θ::Θ
    η::H
    t::T
    a::A
    par :: P
end


@inline function (fam::ExponentialFamily)(par)
    base = fam.base
    θ = fam.θ(par)
    η = fam.η(θ)
    t = fam.t
    a = fam.a(θ)
    ExpFamMeasure(base, θ, η, t, a, par)
end

basemeasure(d::ExpFamMeasure) = d.base

function logdensity_def(d::ExpFamMeasure, x)
    mydot(d.η, d.t(x)) - d.a
end


function MeasureBase.powermeasure(fam::ExponentialFamily, dims::NTuple{N,I}) where {N,I}
    @inline t(x) = sum(fam.t, x)
    a = AffineTransform((σ=prod(dims),)) ∘ fam.a
    ExponentialFamily(fam.base ^ dims, fam.θ, fam.η, t, a)
end

struct ExpFamLikelihood{C,Θ,H,T,A} <: AbstractLikelihood
    c::C
    θ::Θ
    η::H
    t::T
    a::A
end

export likelihood

function likelihood(fam::ExponentialFamily, x)
    c = logdensityof(fam.base, x)
    t = fam.t(x)
    ExpFamLikelihood(c, fam.θ, fam.η, t, fam.a)
end

@inline function logdensity_def(ℓ::ExpFamLikelihood, par)
    θ = ℓ.θ(par)
    mydot(ℓ.η(θ), ℓ.t) - ℓ.a(θ) + ℓ.c 
end

basemeasure(fam::ExponentialFamily) = fam.base