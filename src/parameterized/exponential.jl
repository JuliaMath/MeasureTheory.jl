
# Exponential distribution

import StatsFuns
export Exponential

@parameterized Exponential(λ) ≃ Lebesgue(ℝ₊)

@kwstruct Exponential(λ)
@kwstruct Exponential()

function logdensity(d::Exponential{()} , x)
    return -x
end

Base.rand(rng::AbstractRNG, T::Type, μ::Exponential{()}) = randexp(rng,T)


##########################

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:λ,)})
    randexp(rng, T) / d.λ
end

function logdensity(d::Exponential{(:λ,)}, x)
    z = x * d.λ
    return logdensity(Exponential(), z) + log(d.λ)
end

distproxy(d::Exponential{(:λ,)}) = Dists.Exponential(d.λ)
