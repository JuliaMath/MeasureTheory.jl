Base.rand(μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, μ)

function Base.rand(rng::AbstractRNG, d::ProductMeasure)
    tuple((rand(rng, dn) for dn in d.data)...)
end


@inline Random.rand!(d::AbstractMeasure, arr::AbstractArray) = rand!(GLOBAL_RNG, d, arr)


# struct ArraySlot{A,I}
#     arr::A
#     i::I
# end

# function rand!(rng::AbstractRNG, d::AbstractMeasure, x::ArraySlot)
#     x.arr[x.i...] = rand(rng, d)
# end
