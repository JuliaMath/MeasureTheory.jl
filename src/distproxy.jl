export proxy
function proxy end

import Statistics
import StatsBase

PROXIES = Dict(
    :Statistics => [
        :mean
        :std
        :var
        :quantile
    ],
    :StatsBase => [:entropy],
    :Distributions => [
        :cdf
        :ccdf
        :logcdf
        :logccdf
    ],
)

for m in keys(PROXIES)
    for f in PROXIES[m]
        @eval begin
            import $m: $f
            export $f
        end
    end
end

entropy(m::AbstractMeasure, b::Real) = entropy(proxy(m), b)
mean(m::AbstractMeasure) = mean(proxy(m))
std(m::AbstractMeasure) = std(proxy(m))
var(m::AbstractMeasure) = var(proxy(m))
quantile(m::AbstractMeasure, q) = quantile(proxy(m), q)

for f in [
    :cdf
    :ccdf
    :logcdf
    :logccdf
]
    @eval begin
        $f(d::AbstractMeasure, args...) = $f(MeasureTheory.proxy(d), args...)
    end
end
