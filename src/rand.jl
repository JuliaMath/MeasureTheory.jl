Base.rand(μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, sampletype(μ), μ)

Base.rand(T::Type, d::AbstractMeasure) = rand(GLOBAL_RNG, T, d)

Base.rand(rng::AbstractRNG, d::AbstractMeasure) = rand(rng, sampletype(d), d)

function Base.rand(rng::AbstractRNG, T::Type{<:Tuple}, d::ProductMeasure)
    s = schema(T)
    tuple((rand(rng, sn, dn) for (sn,dn) in zip(s,d.data))...)
end

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
