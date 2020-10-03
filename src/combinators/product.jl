export ProductMeasure
using Base: eltype

struct ProductMeasure{T} <: AbstractMeasure
    components :: T

    ProductMeasure(μs...) = new{Tuple{sampletype.(μs)...}}(μs)
end


# ProductMeasure(m::NTuple{N, Measure{X}}) where {N,X} = ProductMeasure(m...)

Base.length(m::ProductMeasure{T}) where {T} = length(m.components)

function Base.:*(μ::ProductMeasure{X}, ν::ProductMeasure{Y}) where {X,Y}
    components = (μ.components..., ν.components...)
    ProductMeasure(components...)
end

function Base.:*(μ, ν::ProductMeasure{Y}) where {Y}
    components = (μ, ν.components...)
    ProductMeasure(components...)
end

function Base.:*(μ::ProductMeasure{X}, ν) where {X}
    components = (μ.components..., ν)
    ProductMeasure(components...)
end

function Base.:*(μ, ν)
    components = (μ, ν)
    ProductMeasure(components...)
end

function Base.rand(μ::ProductMeasure{T}) where T
    return rand.(μ.components)
end

sampletype(μ::ProductMeasure) = Tuple{sampletype.(μ.components)...}
