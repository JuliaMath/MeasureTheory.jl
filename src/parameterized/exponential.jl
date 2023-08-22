
# Exponential distribution

export Exponential

@parameterized Exponential(β)

insupport(::Exponential, x) = x ≥ 0
basemeasure(::Exponential) = LebesgueBase()

@kwstruct Exponential()

@inline function logdensity_def(d::Exponential{()}, x)
    return -x
end

Base.rand(rng::AbstractRNG, T::Type, μ::Exponential{()}) = randexp(rng, T)

##########################
# Scale β

@kwstruct Exponential(β)

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:β,)})
    randexp(rng, T) * d.β
end

@inline function logdensity_def(d::Exponential{(:β,)}, x)
    z = x / d.β
    return logdensity_def(Exponential(), z) - log(d.β)
end

proxy(d::Exponential{(:β,)}) = Dists.Exponential(d.β)

##########################
# Log-Scale logβ

@kwstruct Exponential(logβ)

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:logβ,)})
    randexp(rng, T) * exp(d.logβ)
end

@inline function logdensity_def(d::Exponential{(:logβ,)}, x)
    z = x * exp(-d.logβ)
    return logdensity_def(Exponential(), z) - d.logβ
end

proxy(d::Exponential{(:logβ,)}) = Dists.Exponential(exp(d.logβ))

##########################
# Rate λ

@kwstruct Exponential(λ)

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:λ,)})
    randexp(rng, T) / d.λ
end

@inline function logdensity_def(d::Exponential{(:λ,)}, x)
    z = x * d.λ
    return logdensity_def(Exponential(), z) + log(d.λ)
end

proxy(d::Exponential{(:λ,)}) = Dists.Exponential(1 / d.λ)

##########################
# Log-Rate logλ

@kwstruct Exponential(logλ)

function Base.rand(rng::AbstractRNG, T::Type, d::Exponential{(:logλ,)})
    randexp(rng, T) * exp(-d.logλ)
end

@inline function logdensity_def(d::Exponential{(:logλ,)}, x)
    z = x * exp(d.logλ)
    return logdensity_def(Exponential(), z) + d.logλ
end

proxy(d::Exponential{(:logλ,)}) = Dists.Exponential(exp(-d.logλ))

smf(::Exponential{()}, x) = smf(StdExponential(), x)

invsmf(::Exponential{()}, p) = invsmf(StdExponential(), p)
