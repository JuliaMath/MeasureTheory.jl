export PowerMeasure
using Base: eltype

struct PowerMeasure{M,N}
    μ::M
    size::NTuple{N,Int}
end
ℝⁿ
function Base.:^(μ::Measure{X}, n::Integer) where {X}
    components = ntuple(i -> μ, n)
    ProductMeasure(components...)
end
