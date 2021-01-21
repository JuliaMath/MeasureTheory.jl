Base.rand(μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, sampletype(μ), μ)

# struct ArraySlot{A,I}
#     arr::A
#     i::I
# end

# function rand!(rng::AbstractRNG, d::AbstractMeasure, x::ArraySlot)
#     x.arr[x.i...] = rand(rng, d)
# end
