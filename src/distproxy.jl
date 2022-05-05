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
    :StatsBase => [
        :entropy
    ],
    :Distributions => [
        :cdf
        :ccdf
        :logcdf
        :logccdf
    ]
)

for m in keys(PROXIES)
    for f in PROXIES[m]
        @eval begin
            import $m: $f
            export $f
            $m.$f(d::AbstractMeasure, args...) = $m.$f(MeasureTheory.proxy(d), args...)
        end
    end
end
