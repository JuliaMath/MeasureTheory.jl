struct Truncated{P,L,U,W} <: AbstractMeasure
    base::P
    lower::L
    upper::U
    logweight::W

    function Truncated(d::P, lower::L, upper::U) where {P,L,U}
        tails = cdf(d, lower) + ccdf(d, upper)
        logweight = -log1p(-tails)
        new{P,L,U,typeof(logweight)}(d, lower, upper, logweight)
    end

    function Truncated(d::P, ::Nothing, upper::U) where {P,U}
        logweight = -logcdf(d, upper)
        new{P,Nothing,U,typeof(logweight)}(d, nothing, upper, logweight)
    end

    function Truncated(d::P, lower::L, ::Nothing) where {P,L}
        logweight = -logccdf(d, lower)
        new{P,L,Nothing,typeof(logweight)}(d, lower, nothing, logweight)
    end
end

insupport(d::Truncated, x) = insupport(d.base, x) && d.lower < x < d.upper

insupport(d::Truncated{P,L,Nothing}, x) where {P,L} = insupport(d.base, x) && d.lower < x

insupport(d::Truncated{P,Nothing,U}, x) where {P,U} = insupport(d.base, x) && x < d.upper

export truncated

truncated(d, lower, upper) = Truncated(d, lower, upper)

truncated(d, ::Nothing, ::Nothing) = d

truncated(d; lower=nothing, upper=nothing) = truncated(d, lower, upper)



logdensity_def(d::Truncated, x) = logdensity_def(d.base, x)

basemeasure(d::Truncated, x) = weightedmeasure(d.logweight, basemeasure(d.base, x))
