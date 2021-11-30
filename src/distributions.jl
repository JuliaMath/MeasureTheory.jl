import Distributions

export Dists
const Dists = Distributions

function Measure(dist::Dists.Distribution{F,S}) where {F,S}
    X = Base.eltype(dist)
    DistributionMeasure{F,S,X}(dist)
end

struct DistributionMeasure{F,S,X} <: AbstractMeasure
    dist::Dists.Distribution{F,S}
end

@inline function basemeasure(μ::DistributionMeasure{F,S,X}) where {F,S<:Dists.Continuous,X}
    return Lebesgue(X)
end

export logdensity

@inline function logdensity(μ::DistributionMeasure{F,S,X}, x::X) where {F,S,X}
    return Dists.logpdf(μ.dist, x)
end

function Base.rand(rng::AbstractRNG, μ::DistributionMeasure{F,S,X}) where {F,S,X}
    return rand(rng, μ.dist)
end

basemeasure(d::Dists.AbstractMvNormal) = Lebesgue(ℝ)^size(d)

@inline function basemeasure(μ::Dists.Distribution{Dists.Univariate,Dists.Continuous})
    return Lebesgue(ℝ)
end

@inline function basemeasure(μ::Dists.Distribution{Dists.Univariate,Dists.Discrete})
    return CountingMeasure(ℤ)
end

∫(::typeof(identity), ::Dists.Distribution) = 1.0

logdensity(μ::Dists.Distribution, x) = Dists.logpdf(μ, x)

density(μ::Dists.Distribution, x) = Dists.pdf(μ, x)

testvalue(d::Dists.Distribution) = rand(d)
