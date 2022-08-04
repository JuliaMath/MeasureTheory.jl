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

@inline function logdensity_def(μ::DistributionMeasure{F,S,X}, x::X) where {F,S,X}
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

logdensity_def(μ::Dists.Distribution, x) = Dists.logpdf(μ, x)

density_def(μ::Dists.Distribution, x) = Dists.pdf(μ, x)

testvalue(d::Dists.Distribution) = rand(d)

insupport(d::Dists.Distribution, x) = logdensityof(d, x) > -Inf



@inline function as(d, _data::NamedTuple)
    if hasmethod(Dists.support, (typeof(d),))
        return asTransform(Dists.support(d)) 
    end

    error("Not implemented:\nas($d)")
end

using TransformVariables: ShiftedExp, ScaledShiftedLogistic

function asTransform(supp:: Dists.RealInterval) 
    (lb, ub) = (supp.lb, supp.ub)

    (lb, ub) == (-Inf, Inf) && (return asℝ)
    isinf(ub) && return ShiftedExp(true,lb)
    isinf(lb) && return ShiftedExp(false,lb)
    return ScaledShiftedLogistic(ub-lb, lb)
end

as(μ::AbstractMeasure,  _data::NamedTuple) = as(μ)

as(d::Dists.AbstractMvNormal, _data::NamedTuple = NamedTuple()) = TV.as(Array, size(d))


function as(d::Dists.Distribution{Dists.Univariate}, _data::NamedTuple=NamedTuple())
    sup = Dists.support(d)
    lo = isinf(sup.lb) ? -TV.∞ : sup.lb
    hi = isinf(sup.ub) ? TV.∞ : sup.ub
    as(Real, lo,hi)
end

function as(d::Dists.Product, _data::NamedTuple=NamedTuple())
    n = length(d)
    v = d.v
    as(Vector, as(v[1]), n)
end
