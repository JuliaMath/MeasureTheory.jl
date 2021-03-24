Base.rand(μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, Float64, μ)

Base.rand(T::Type, d::AbstractMeasure) = rand(GLOBAL_RNG, T, d)

Base.rand(rng::AbstractRNG, d::AbstractMeasure) = rand(rng, sampletype(d), d)

function Base.rand(rng::AbstractRNG, T::Type, d::ProductMeasure)
    s = schema(T)
    tuple((rand(rng, sn, dn) for (sn,dn) in zip(s,d.data))...)
end


@inline Random.rand!(d::AbstractMeasure, arr::AbstractArray) = rand!(GLOBAL_RNG, d, arr)


# struct ArraySlot{A,I}
#     arr::A
#     i::I
# end

# function rand!(rng::AbstractRNG, d::AbstractMeasure, x::ArraySlot)
#     x.arr[x.i...] = rand(rng, d)
# end
