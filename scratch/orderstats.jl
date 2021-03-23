struct OrderStatistic{n,r,B} <: AbstractMeasure
    base::B
end

function basemeasure(d::OrderStatistic{n,r})
    WeightedMeasure(logbeta(n+1, r+1), d.base)
end

function logdensity(d::OrderStatistic{n,r}, x)
    logF(x) = logcdf(d.base, x)
    logFᶜ(x) = logccdf(d.base, x)
    return (r - 1) * logF(x) + (n - r) * logFᶜ(x)
end
