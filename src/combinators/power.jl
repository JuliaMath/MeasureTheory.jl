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
struct PowerMeasure{M,N}
    μ::M
    size::NTuple{N,Int}
end

function Base.:^(μ::Measure{X}, n::Integer) where {X}
    components = ntuple(i -> μ, n)
    ProductMeasure(components...)
end
