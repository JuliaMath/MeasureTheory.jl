
# Exponential distribution

import StatsFuns
export Exponential

@measure Exponential(λ) ≪ Lebesgue(ℝ₊)

function logdensity(d::Exponential{()} , x)
    return -x
end

Base.rand(rng::AbstractRNG, T::Type, μ::Exponential{()}) = -log(rand(rng, T))

≪(::Exponential, ::Lebesgue{ℝ₊}) = true

representative(::Exponential) = Lebesgue(ℝ₊)


##########################

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:min, :λ)})
    rand(rng, T, Exponential()) / d.λ + d.min
end

function logdensity(d::Exponential{(:min, :λ)}, x)
    z = (x - d.min) * d.λ
    return logdensity(Exponential(), z) + log(d.λ)
end

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:λ,)})
    rand(rng, T, Exponential()) / d.λ
end

function logdensity(d::Exponential{(:λ,)}, x)
    z = x * d.λ
    return logdensity(Exponential(), z) + log(d.λ)
end

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:min,)})
    rand(rng, T, Exponential()) + d.min
end

function logdensity(d::Exponential{(:min,)}, x)
    z = x - d.min
    return logdensity(Exponential(), z)
end
