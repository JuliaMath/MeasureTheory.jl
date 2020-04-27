import Distributions

export Dists
const Dists = Distributions

function Measure(dist::Dists.Distribution{F,S}) where {F,S}
    X = eltype(dist)
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
