struct ExponentialFamily{M,H,T,A} <: AbstractKleisli
    parent::M
    η::H
    t::T
    a::A
end

function (fam::ExponentialFamily)(θ)
    η = fam.η(θ)
    a = fam.a(θ)
    ℓ(x) = dot(η, fam.t(x)) - a
    ∫exp(ℓ, fam.parent)
end

struct ExpFamLikelihood{C,H,T,A} <: AbstractLikelihood
    c::C
    η::H
    t::T
    a::A
end


function likelihood(fam::ExponentialFamily, x)
    c = logdensity_def(fam.parent, x)
    t = fam.t(x)
    ExpFamLikelihood(c, fam.η, t, fam.a)
end

function logdensity_def(ℓ::ExpFamLikelihood, θ)
    dot(ℓ.η(θ), ℓ.t) - ℓ.a(θ) + ℓ.c 
end

