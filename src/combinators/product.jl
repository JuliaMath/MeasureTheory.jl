export ProductMeasure

using MappedArrays
using Base: @propagate_inbounds

struct ProductMeasure{T} <: AbstractMeasure
    data :: T

    ProductMeasure(μs...) = new{typeof(μs)}(μs)
    ProductMeasure(μs) = new{typeof(μs)}(μs)
end

Base.size(μ::ProductMeasure) = size(μ.data)

function Base.show(io::IO, μ::ProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, join(string.(μ.data), " ⊗ "))
end

function Base.show_unquoted(io::IO, μ::ProductMeasure, indent::Int, prec::Int)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end

# ProductMeasure(m::NTuple{N, Measure{X}}) where {N,X} = ProductMeasure(m...)

Base.length(m::ProductMeasure{T}) where {T} = length(m.data)

function Base.:*(μ::ProductMeasure{Tuple{}}, ν::N) where {X, N <: AbstractMeasure}
    ProductMeasure((ν,))
end

function Base.:*(μ::ProductMeasure{X}, ν::ProductMeasure{Y}) where {X,Y}
    data = (μ.data..., ν.data...)
    ProductMeasure(data...)
end

function Base.:*(μ, ν::ProductMeasure{Y}) where {Y}
    data = (μ, ν.data...)
    ProductMeasure(data...)
end

function Base.:*(μ::ProductMeasure{X}, ν::N) where {X, N <: AbstractMeasure}
    data = (μ.data..., ν)
    ProductMeasure(data...)
end

function Base.:*(μ::M, ν::N) where {M <: AbstractMeasure, N <: AbstractMeasure}
    data = (μ, ν)
    ProductMeasure(data...)
end

using Tullio

function logdensity(d::ProductMeasure, x::NTuple)
    mapreduce(logdensity, +, d.data, x)
end

function logdensity(d::ProductMeasure{A}, x::NTuple) where {A <: AbstractArray}
    mapreduce(logdensity, +, d.data, x)
end

function logdensity(d::ProductMeasure{A}, x) where {T, A<:AbstractArray}
    data = d.data
    @boundscheck size(data) == size(x) || throw(BoundsError)
    s = 0.0
    Δs(i) = logdensity(data[i], x[i])
    for i in eachindex(x)
        s += Δs(i)
    end
    return s
end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,1}}
#     data = d.data
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = logdensity(data[i], x[i])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,2}}
#     data = d.data
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j], x[i,j])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,3}}
#     data = d.data
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j,k], x[i,j,k])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,4}}
#     data = d.data
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j,k,l], x[i,j,k,l])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,5}}
#     data = d.data
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j,k,l,m], x[i,j,k,l,m])
#     s
# end

export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG

@propagate_inbounds function Random.rand!(rng::AbstractRNG, d::ProductMeasure, x::AbstractArray)
    @boundscheck size(d.data) == size(x) || throw(BoundsError)

    @inbounds for j in eachindex(x)
        x[j] = rand(rng, eltype(x), d.data[j])
    end
    x
end


# function Base.rand(rng::AbstractRNG, d::ProductMeasure)
#     return rand(rng, sampletype(d), d)
# end

# function Base.rand(T::Type, d::ProductMeasure)
#     return rand(Random.GLOBAL_RNG, T, d)
# end

# function Base.rand(d::ProductMeasure)
#     T = sampletype(d)
#     return rand(Random.GLOBAL_RNG, T, d)
# end

function sampletype(d::ProductMeasure{A}) where {T,N,A <: AbstractArray{T,N}}
    S = @inbounds sampletype(d.data[1])
    Array{S, N}
end

function sampletype(d::ProductMeasure{<: Tuple}) 
    Tuple{sampletype.(d.data)...}
end

# TODO: Pull weights outside
basemeasure(μ::ProductMeasure) = ProductMeasure(basemeasure.(μ.data))

# function logdensity(μ::ProductMeasure{Aμ}, x::Ax) where {Aμ <: MappedArray, Ax <: AbstractArray}
#     μ.data
# end
