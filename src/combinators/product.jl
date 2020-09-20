export ProductMeasure
using Base: eltype

struct ProductMeasure{T} <: AbstractMeasure{T}
    # TODO: Type annotation
    components

    ProductMeasure(μs...) = new{Tuple{eltype.(μs)...}}(μs)
end


# ProductMeasure(m::NTuple{N, Measure{X}}) where {N,X} = ProductMeasure(m...)

Base.length(m::ProductMeasure{T}) where {T} = length(m.components)

function Base.:*(μ::ProductMeasure{X}, ν::ProductMeasure{Y}) where {X,Y}
    components = (μ.components..., ν.components...)
    ProductMeasure(components...)
end

function Base.:*(μ::Measure{X}, ν::ProductMeasure{Y}) where {X,Y}
    components = (μ, ν.components...)
    ProductMeasure(components...)
end

function Base.:*(μ::ProductMeasure{X}, ν::Measure{Y}) where {X,Y}
    components = (μ.components..., ν)
    ProductMeasure(components...)
end

function Base.:*(μ::Measure{X}, ν::Measure{Y}) where {X,Y}
    components = (μ, ν)
    ProductMeasure(components...)
end

export × 
function ×(μ::Measure{X}, ν::Measure{Y}) where {X,Y}
    return μ*ν
end

function Base.rand(μ::ProductMeasure{T}) where T
    return rand.(μ.components)
end
