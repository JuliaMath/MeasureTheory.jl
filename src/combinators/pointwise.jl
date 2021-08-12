export ⊙

@concrete terse struct PointwiseProductMeasure{T} <: AbstractMeasure
    data :: T

    PointwiseProductMeasure(μs...) = new{typeof(μs)}(μs)
    PointwiseProductMeasure(μs) = new{typeof(μs)}(μs)
end

Base.size(μ::PointwiseProductMeasure) = size(μ.data)

function Base.show(io::IO, μ::PointwiseProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, join(string.(μ.data), " ⊙ "))
end

function Base.show_unquoted(io::IO, μ::PointwiseProductMeasure, indent::Int, prec::Int)
    io = IOContext(io, :compact => true)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end

Base.length(m::PointwiseProductMeasure{T}) where {T} = length(m.data)

function ⊙(μ::PointwiseProductMeasure{X}, ν::PointwiseProductMeasure{Y}) where {X,Y}
    data = (μ.data..., ν.data...)
    PointwiseProductMeasure(data...)
end

function ⊙(μ::AbstractMeasure, ν::PointwiseProductMeasure)
    data = (μ, ν.data...)
    PointwiseProductMeasure(data...)
end

function ⊙(μ::PointwiseProductMeasure, ν::N) where {N <: AbstractMeasure}
    data = (μ.data..., ν)
    PointwiseProductMeasure(data...)
end

function ⊙(μ::M, ν::N) where {M <: AbstractMeasure, N <: AbstractMeasure}
    data = (μ, ν)
    PointwiseProductMeasure(data...)
end

function ⊙(μ::AbstractMeasure, ℓ::Likelihood)
    data = (μ, ℓ)
    PointwiseProductMeasure(data...)
end

function logdensity(d::PointwiseProductMeasure, x)
    sum((logdensity(dⱼ, x) for dⱼ in d.data))
end

function sampletype(d::PointwiseProductMeasure) 
    @inbounds sampletype(first(d.data))
end

basemeasure(d::PointwiseProductMeasure) =  @inbounds basemeasure(first(d.data))
