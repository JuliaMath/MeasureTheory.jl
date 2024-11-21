struct Censored{P,L,U,W} <: AbstractMeasure
    parent::P
    lower::L
    upper::U
    ℓ_lower
    ℓ_upper

    function Censored(d::P, lower::L, upper::U) where {P,L,U}
        tails = cdf(d, lower) + ccdf(d, upper)
        logweight = -log1p(-tails)
        new{P,L,U,typeof(logweight)}(d, lower, upper, logweight)
    end

    function Censored(d::P, ::Nothing, upper::U) where {P,U}
        logweight = -logcdf(d, upper)
        new{P,Nothing,U,typeof(logweight)}(d, nothing, upper, logweight)
    end

    function Censored(d::P, lower::L, ::Nothing) where {P,L}
        logweight = -logccdf(d, lower)
        new{P,L,Nothing,typeof(logweight)}(d, lower, nothing, logweight)
    end
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::Censored) where {T}
    x = rand(rng, T, d.parent)
    clamp(x, d.lower, d.upper)
end

insupport(d::Censored, x) = insupport(d.parent, x) && d.lower ≤ x ≤ d.upper

insupport(d::Censored{P,L,Nothing}, x) where {P,L} = insupport(d.parent, x) && d.lower ≤ x

insupport(d::Censored{P,Nothing,U}, x) where {P,U} = insupport(d.parent, x) && x ≤ d.upper

export censored

censored(d, lower, upper) = Censored(d, lower, upper)

censored(d, ::Nothing, ::Nothing) = d

censored(d; lower=nothing, upper=nothing) = censored(d, lower, upper)

logdensity_def(d::Censored, x) = logdensity_def(d.parent, x)

basemeasure(d::Censored, x) = basemeasure(d.parent)
