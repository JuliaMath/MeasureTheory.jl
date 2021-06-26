Base.rand(d::AbstractMeasure) = rand(Random.GLOBAL_RNG, Float64, d)

Base.rand(T::Type, μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, T, μ)

Base.rand(rng::AbstractRNG, d::AbstractMeasure) = rand(rng, Float64, d)

Base.rand(rng::AbstractRNG, T::Type, d::ParameterizedMeasure) = rand(rng, distproxy(d))

@inline Random.rand!(d::AbstractMeasure, arr::AbstractArray) = rand!(GLOBAL_RNG, d, arr)



# struct ArraySlot{A,I}
#     arr::A
#     i::I
# end

# function rand!(rng::AbstractRNG, d::AbstractMeasure, x::ArraySlot)
#     x.arr[x.i...] = rand(rng, d)
# end
