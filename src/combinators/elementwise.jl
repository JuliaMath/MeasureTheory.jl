export ⊙

struct ElementwiseProductMeasure{T} <: AbstractMeasure
    data :: T

    ElementwiseProductMeasure(μs...) = new{typeof(μs)}(μs)
    ElementwiseProductMeasure(μs) = new{typeof(μs)}(μs)
end

Base.size(μ::ElementwiseProductMeasure) = size(μ.data)

function Base.show(io::IO, μ::ElementwiseProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, join(string.(μ.data), " ⊙ "))
end

function Base.show_unquoted(io::IO, μ::ElementwiseProductMeasure, indent::Int, prec::Int)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end

Base.length(m::ElementwiseProductMeasure{T}) where {T} = length(m.data)

function ⊙(μ::ElementwiseProductMeasure{X}, ν::ElementwiseProductMeasure{Y}) where {X,Y}
    data = (μ.data..., ν.data...)
    ElementwiseProductMeasure(data...)
end

function ⊙(μ, ν::ElementwiseProductMeasure{Y}) where {Y}
    data = (μ, ν.data...)
    ElementwiseProductMeasure(data...)
end

function ⊙(μ::ElementwiseProductMeasure{X}, ν::N) where {X, N <: AbstractMeasure}
    data = (μ.data..., ν)
    ElementwiseProductMeasure(data...)
end

function ⊙(μ::M, ν::N) where {M <: AbstractMeasure, N <: AbstractMeasure}
    data = (μ, ν)
    ElementwiseProductMeasure(data...)
end

function logdensity(d::ElementwiseProductMeasure, x)
    sum((logdensity(dⱼ, x) for dⱼ in d.data))
end

function sampletype(d::ElementwiseProductMeasure) 
    @inbounds sampletype(first(d.data))
end

basemeasure(μ::ElementwiseProductMeasure) =  @inbounds basemeasure(first(d.data))
