export ProductMeasure

using MappedArrays
using Base: @propagate_inbounds

struct ProductMeasure{F,I} <: AbstractMeasure
    f::F
    inds::I
end

Base.size(μ::ProductMeasure) = size(marginals(μ))

Base.length(m::ProductMeasure{T}) where {T} = length(marginals(μ))

# TODO: Pull weights outside
basemeasure(d::ProductMeasure) = ProductMeasure(basemeasure ∘ d.f, d.inds)


export marginals


testvalue(d::ProductMeasure) = map(testvalue, marginals(d))


###############################################################################
# I <: Tuple

export ⊗
⊗(μs::AbstractMeasure...) = ProductMeasure(identity, μs)

marginals(d::ProductMeasure{F,T}) where {F, T<:Tuple} = map(d.f, d.inds)

function Base.show(io::IO, μ::ProductMeasure{F,T}) where {F,T <: Tuple}
    io = IOContext(io, :compact => true)
    print(io, join(string.(marginals(μ)), " ⊗ "))
end

function logdensity(d::ProductMeasure{F,T}, x::Tuple) where {F,T<:Tuple}
    mapreduce(logdensity, +, d.f.(d.inds), x)
end

###############################################################################
# I <: AbstractArray

marginals(d::ProductMeasure{F,A}) where {F,A<:AbstractArray} = mappedarray(d.f, d.inds)


function logdensity(d::ProductMeasure, x::AbstractArray)
    mar = marginals(d)
    @boundscheck size(mar) == size(x) || throw(BoundsError)
    s = 0.0
    Δs(i) = logdensity(mar[i], x[i])
    for i in eachindex(x)
        s += Δs(i)
    end
    return s
end


###############################################################################
# I <: Base.Generator

###############################################################################

function Base.show(io::IO, μ::ProductMeasure{NamedTuple{N,T}}) where {N,T}
    io = IOContext(io, :compact => true)
    print(io, "Product(",μ.data, ")")
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





using Tullio




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


# function logdensity(μ::ProductMeasure{Aμ}, x::Ax) where {Aμ <: MappedArray, Ax <: AbstractArray}
#     μ.data
# end
