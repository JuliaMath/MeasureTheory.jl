export ProductMeasure

using MappedArrays
using FillArrays
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

@propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure, x)
    data = d.data
    @boundscheck size(data) == size(x) || throw(BoundsError)

    s = 0.0
    Δs(j) = @inbounds logdensity(data[j], x[j])

    @inbounds @simd for j in eachindex(x)
        s += Δs(j)
    end
    s
end

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

# TODO: Pull weights outside
basemeasure(μ::ProductMeasure) = ProductMeasure(basemeasure.(μ.data))

# function logdensity(μ::ProductMeasure{Aμ}, x::Ax) where {Aμ <: MappedArray, Ax <: AbstractArray}
#     μ.data
# end
