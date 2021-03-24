
# Exponential distribution

import StatsFuns
export Exponential

@measure Exponential(λ)

basemeasure(::Exponential) = Lebesgue(ℝ₊)

function logdensity(d::Exponential{()} , x)
    return -x
end

Base.rand(rng::AbstractRNG, μ::Exponential{()}) = randexp(rng)

representative(::Exponential) = Lebesgue(ℝ₊)


##########################

function Base.rand(rng::AbstractRNG, d::Exponential{(:λ,)})
    randexp(rng) / d.λ
end

function logdensity(d::Exponential{(:λ,)}, x)
    z = x * d.λ
    return logdensity(Exponential(), z) + log(d.λ)
end

distproxy(d::Exponential{(:λ,)}) = Dists.Exponential(d.λ)
