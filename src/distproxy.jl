export proxy
function proxy end

PROXIES = Dict(:Distributions => [
    :mean
    :std
    :entropy
    :cdf
])

for m in keys(PROXIES)
    for f in PROXIES[m]
        @eval begin
            import $m: $f
            export $f
            $m.$f(d::AbstractMeasure, args...) = $m.$f(MeasureTheory.proxy(d), args...)
        end
    end
end
