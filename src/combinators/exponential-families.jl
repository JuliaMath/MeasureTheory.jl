export ExponentialFamily

struct ExponentialFamily{B,Θ,H,T,A} <: AbstractKleisli
    base::B
    θ::Θ
    η::H
    t::T
    a::A
end


function MeasureBase.powermeasure(fam::ExponentialFamily, dims::NTuple{N,I}) where {N,I}
    @inline t(x) = sum(fam.t, x)
    a = AffineTransform((σ=prod(dims),)) ∘ fam.a
    ExponentialFamily(fam.base ^ dims, fam.θ, fam.η, t, a)
end

function (fam::ExponentialFamily)(par)
    θ = fam.θ(par)
    η = fam.η(θ)
    a = fam.a(θ)
    @inline ℓ(x) = mydot(η, fam.t(x)) - a
    ∫exp(ℓ, fam.base)
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
