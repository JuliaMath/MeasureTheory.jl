export PowerMeasure
using Base: eltype

"""
    struct PowerMeasure{M,N}
        μ::M
        size::NTuple{N,Int}
    end
    
A power measure is a product of a measure with itself. The number of elements in
the product determines the dimensionality of the resulting support.

Note that power measures are only well-defined for integer powers.

The nth power of a measure μ can be written μ^x.
"""
struct PowerMeasure{M,N} <: AbstractMeasure
    μ::M
    size::NTuple{N,Int}
end

import Base

Base.show(io::IO, μ::PowerMeasure) = print(io, μ.μ, " ^ ", μ.size)

function Base.:^(μ::AbstractMeasure, n::Integer)  
    PowerMeasure(μ, (n,))
end

function Base.:^(μ::AbstractMeasure, size::NTuple{N,I}) where {N, I <: Integer} 
    PowerMeasure(μ, size)
end

sampletype(d::PowerMeasure{M,N}) where {M,N} = AbstractArray{sampletype(d.μ), N}

export rand!

function rand!(result::AbstractArray, d::PowerMeasure)
    @inbounds for j in eachindex(result)
        result[j] = rand(d.μ)
    end
    return result
end

function Base.rand(d::PowerMeasure)
    result = Array{sampletype(d.μ), length(d.size)}(undef, d.size...)
    return rand!(result, d)
end    

basemeasure(μ::PowerMeasure) = basemeasure(μ.μ)^μ.size
