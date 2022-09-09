import Distributions

export Dists
const Dists = Distributions

@inline function as(d, _data::NamedTuple)
    if hasmethod(Dists.support, (typeof(d),))
        return asTransform(Dists.support(d))
    end

    error("Not implemented:\nas($d)")
end

using TransformVariables: ShiftedExp, ScaledShiftedLogistic

function asTransform(supp::Dists.RealInterval)
    (lb, ub) = (supp.lb, supp.ub)

    (lb, ub) == (-Inf, Inf) && (return asℝ)
    isinf(ub) && return ShiftedExp(true, lb)
    isinf(lb) && return ShiftedExp(false, lb)
    return ScaledShiftedLogistic(ub - lb, lb)
end

as(μ::AbstractMeasure, _data::NamedTuple) = as(μ)

as(d::Dists.AbstractMvNormal, _data::NamedTuple = NamedTuple()) = TV.as(Array, size(d))

function as(d::Dists.Distribution{Dists.Univariate}, _data::NamedTuple = NamedTuple())
    sup = Dists.support(d)
    lo = isinf(sup.lb) ? -TV.∞ : sup.lb
    hi = isinf(sup.ub) ? TV.∞ : sup.ub
    as(Real, lo, hi)
end

function as(d::Dists.Product, _data::NamedTuple = NamedTuple())
    n = length(d)
    v = d.v
    as(Vector, as(v[1]), n)
end

as(m::DistributionMeasures.DistributionMeasure) = as(m.d)
