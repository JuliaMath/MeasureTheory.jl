import Distributions

export Dists
const Dists = Distributions


function Measure(dist::Dists.Distribution{F,S}) where {F,S}
    X = Base.eltype(dist)
    DistributionMeasure{F,S,X}(dist)
end

struct DistributionMeasure{F, S, X} <: AbstractMeasure
    dist :: Dists.Distribution{F, S}
end

function basemeasure(μ::DistributionMeasure{F,S,X}) where {F, S <: Dists.Continuous, X}
    return Lebesgue(X)
end

export logdensity

function logdensity(μ::DistributionMeasure{F,S,X}, x::X) where {F, S, X}
    return Dists.logpdf(μ.dist, x)
end

function Base.rand(μ::DistributionMeasure{F,S,X}) where {F, S, X}
    return rand(μ.dist)
end

export Normal




function basemeasure(μ::Dists.Distribution{Dists.Univariate,Dists.Continuous})
    T = eltype(μ)
    return Lebesgue(T)
end

function basemeasure(μ::Dists.Distribution{Dists.Multivariate,Dists.Continuous})
    x = rand(μ)
    return Lebesgue(typeof(x))
end

logdensity(μ::Dists.Distribution, x) = Dists.logpdf(μ,x)

density(μ::Dists.Distribution, x) = Dists.pdf(μ,x)
