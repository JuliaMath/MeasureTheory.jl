import Distributions

export Dists
const Dists = Distributions

function Measure(dist::Dists.Distribution{F,S}) where {F,S}
    X = Base.eltype(dist)
    DistributionMeasure{F,S,X}(dist)
end

struct DistributionMeasure{F, S, X} <: Measure{X}
    dist :: Dists.Distribution{F, S}
end

function baseMeasure(μ::DistributionMeasure{F,S,X}) where {F, S <: Dists.Continuous, X}
    return Lebesgue(X)
end

export logdensity

function logdensity(μ::DistributionMeasure{F,S,X}, x::X) where {F, S, X}
    return Dists.logpdf(μ.dist, x)
end

export rand

function Base.rand(μ::DistributionMeasure{F,S,X}) where {F, S, X}
    return rand(μ.dist)
end

export Normal

function Normal(μ::T, σ::T) where {T <: AbstractFloat}
    d = Dists.Normal(μ,σ; check_args=false)
    return Measure(d)
end

Normal(μ::Real, σ::Real) = Normal(promote(μ, σ)...)
Normal(μ::Integer, σ::Integer) = Normal(float(μ), float(σ))
Normal(μ::T) where {T <: Real} = Normal(μ, one(T))
