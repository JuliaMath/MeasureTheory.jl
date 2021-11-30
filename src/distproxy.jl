export distproxy
function distproxy end

# import MonteCarloMeasurements

PROXIES = Dict(
    :Distributions => [
        :mean
        :std
        :entropy
        :cdf
    ],
    # :MonteCarloMeasurements => [
    #     :Particles
    # ]
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

Base.rand(rng::AbstractRNG, ::Type{T}, d::ParameterizedMeasure) where {T} =
    rand(rng, T, proxy(rand, d))

# MonteCarloMeasurements.Particles(N::Int, d::AbstractMeasure) = MonteCarloMeasurements.Particles(N, distproxy(d))

# using MonteCaroMeasurements

# MonteCaroMeasurementsPROXIES = [
#     :Particles
# ]

# for f in DistributionsPROXIES
#     @eval begin
#         import Distributions: $f
#         export $f
#         Distributions.$f(d::AbstractMeasure) = Distributions.$f(MeasureTheory.distproxy(d))
#     end
# end
