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
    @boundscheck size(d.inds) == size(x) || throw(BoundsError)
    mar = marginals(d)
    s = 0.0
    Δs(i) = logdensity(mar[i], x[i])
    for i in eachindex(x)
        s += Δs(i)
    end
    return s
end


function TV.as(d::ProductMeasure{F,A}) where {F,A<:AbstractArray}
    d1 = marginals(d).f(first(marginals(d).data))
    as(Array, as(d1), size(marginals(d))...)
end

function Base.show(io::IO, d::ProductMeasure{F,A}) where {F,A<:AbstractArray}
    print(io, "For(")
    print(io, d.f, ", ")
    print(io, d.inds, ")")
end


###############################################################################
# I <: CartesianIndices

function Base.show(io::IO, d::ProductMeasure{F,I}) where {F, I<:CartesianIndices}
    print(io, "For(")
    print(io, d.f, ", ")
    join(io, size(d.inds), ", ")
    print(io, ")")
end

###############################################################################
# I <: Base.Generator

function TV.as(d::ProductMeasure{F,I}) where {F, I<:Base.Generator}
    d1 = marginals(d).f(first(marginals(d).iter))
    as(Array, as(d1), size(marginals(d))...) 
end



# TODO Make this reproducible
function Base.rand(rng::AbstractRNG, T::Type, d::ProductMeasure{F,I}) where {F,I<:Base.Generator}
    r(x) = rand(rng, T, x)
    Base.Generator(r ∘ marginals(d).f, marginals(d).iter)
end

function logdensity(d::ProductMeasure{F,I}, x) where {F, I<:Base.Generator}
    sum((logdensity(dj, xj) for (dj, xj) in zip(marginals(d), x)))
end






# ProductMeasure(m::NTuple{N, Measure{X}}) where {N,X} = ProductMeasure(m...)





using Tullio




# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,1}}
#     data = marginals(d)
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = logdensity(data[i], x[i])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,2}}
#     data = marginals(d)
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j], x[i,j])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,3}}
#     data = marginals(d)
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j,k], x[i,j,k])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,4}}
#     data = marginals(d)
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j,k,l], x[i,j,k,l])
#     s
# end

# @propagate_inbounds function MeasureTheory.logdensity(d::ProductMeasure{A}, x) where{T, A<:AbstractArray{T,5}}
#     data = marginals(d)
#     @boundscheck size(data) == size(x) || throw(BoundsError)
#     @tullio s = @inbounds logdensity(data[i,j,k,l,m], x[i,j,k,l,m])
#     s
# end

export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG

@propagate_inbounds function Random.rand!(rng::AbstractRNG, d::ProductMeasure, x::AbstractArray)
    @boundscheck size(marginals(d)) == size(x) || throw(BoundsError)

    @inbounds for j in eachindex(x)
        x[j] = rand(rng, eltype(x), marginals(d)[j])
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
    S = @inbounds sampletype(marginals(d)[1])
    Array{S, N}
end

function sampletype(d::ProductMeasure{<: Tuple}) 
    Tuple{sampletype.(marginals(d))...}
end


# function logdensity(μ::ProductMeasure{Aμ}, x::Ax) where {Aμ <: MappedArray, Ax <: AbstractArray}
#     μ.data
# end
