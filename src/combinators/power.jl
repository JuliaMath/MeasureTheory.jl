import Base

# """
# A power measure is a product of a measure with itself. The number of elements in
# the product determines the dimensionality of the resulting support.

# Note that power measures are only well-defined for integer powers.

# The nth power of a measure μ can be written μ^x.
# """
# PowerMeasure{M,N,D} = ProductMeasure{Fill{M,N,D}}

# function Base.show(io::IO, μ::PowerMeasure)
#     io = IOContext(io, :compact => true)
#     print(io, μ.data.value, " ^ ", size(μ.data))
# end

# function Base.show_unquoted(io::IO, μ::PowerMeasure{M,N,D}, indent::Int, prec::Int) where {M,N,D}
#     io = IOContext(io, :compact => true)
#     if Base.operator_precedence(:^) ≤ prec
#         print(io, "(")
#         show(io, μ.data.value)
#         print(io, ")")
#     else
#         show(io, size(μ.data))
#     end
#     return nothing
# end

export PowerMeasure

PowerMeasure{D,I} = ProductMeasure{Const{D},I}

function Base.:^(μ::AbstractMeasure, dims::Integer...) where {N, I<:Integer}
    return μ^dims
end

function Base.:^(μ::M, dims::NTuple{N,I}) where {M <: AbstractMeasure, N, I<:Integer}
    C = constructorof(M)
    p = params(μ)

    For(C, Fill(p, dims...))
end


# sampletype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{sampletype(first(marginals(d))), N}

function Base.:^(μ::WeightedMeasure, dims::NTuple{N,I}) where {N, I<:Integer}
    k = prod(dims) * μ.logweight
    return WeightedMeasure(k, μ.base^dims)
end

params(::PowerMeasure{D})       where {D} = params(D)
params(::Type{PowerMeasure{D}}) where {D} = params(D)
    
# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)
