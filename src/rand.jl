Base.rand(μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, sampletype(μ), μ)

function Base.rand(rng::AbstractRNG, ::Type{A}, d::AbstractMeasure) where {T, N, A <: AbstractArray{T,N}}
    dims = size(d)
    x = A(undef, dims...)
    return @inbounds rand!(rng, d, x)
end

@inline Random.rand!(d::AbstractMeasure, arr::AbstractArray) = rand!(GLOBAL_RNG, d, arr)


# struct ArraySlot{A,I}
#     arr::A
#     i::I
# end

# function rand!(rng::AbstractRNG, d::AbstractMeasure, x::ArraySlot)
#     x.arr[x.i...] = rand(rng, d)
# end
