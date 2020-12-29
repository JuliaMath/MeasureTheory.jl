export ProductMeasure
using Base: eltype

struct ProductMeasure{T} <: AbstractMeasure
    components :: T

    ProductMeasure(μs...) = new{typeof(μs)}(μs)
    ProductMeasure(μs) = new{typeof(μs)}(μs)
end

function Base.show(io::IO, μ::ProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, join(string.(μ.components), " ⊗ "))
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

Base.length(m::ProductMeasure{T}) where {T} = length(m.components)

function Base.:*(μ::ProductMeasure{X}, ν::ProductMeasure{Y}) where {X,Y}
    components = (μ.components..., ν.components...)
    ProductMeasure(components...)
end

function Base.:*(μ, ν::ProductMeasure{Y}) where {Y}
    components = (μ, ν.components...)
    ProductMeasure(components...)
end

function Base.:*(μ::ProductMeasure{X}, ν::N) where {X, N <: AbstractMeasure}
    components = (μ.components..., ν)
    ProductMeasure(components...)
end

function Base.:*(μ::M, ν::N) where {M <: AbstractMeasure, N <: AbstractMeasure}
    components = (μ, ν)
    ProductMeasure(components...)
end

function Base.rand(μ::ProductMeasure{T}) where T
    return rand.(μ.components)
end

sampletype(μ::ProductMeasure) = Tuple{sampletype.(μ.components)...}

basemeasure(μ::ProductMeasure) = ProductMeasure(basemeasure.(μ.components))
