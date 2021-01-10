export ProductMeasure
using Base: eltype

struct ProductMeasure{T} <: AbstractMeasure
    data :: T

    ProductMeasure(μs...) = new{typeof(μs)}(μs)
    ProductMeasure(μs) = new{typeof(μs)}(μs)
end

function Base.show(io::IO, μ::ProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, join(string.(μ.data), " ⊗ "))
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

Base.length(m::ProductMeasure{T}) where {T} = length(m.data)

function Base.:*(μ::ProductMeasure{X}, ν::ProductMeasure{Y}) where {X,Y}
    data = (μ.data..., ν.data...)
    ProductMeasure(data...)
end

function Base.:*(μ, ν::ProductMeasure{Y}) where {Y}
    data = (μ, ν.data...)
    ProductMeasure(data...)
end

function Base.:*(μ::ProductMeasure{X}, ν::N) where {X, N <: AbstractMeasure}
    data = (μ.data..., ν)
    ProductMeasure(data...)
end

function Base.:*(μ::M, ν::N) where {M <: AbstractMeasure, N <: AbstractMeasure}
    data = (μ, ν)
    ProductMeasure(data...)
end

function Base.rand(μ::ProductMeasure{T}) where T
    return rand.(μ.data)
end

sampletype(μ::ProductMeasure) = Tuple{sampletype.(μ.data)...}

basemeasure(μ::ProductMeasure) = ProductMeasure(basemeasure.(μ.data))

# function logdensity(μ::ProductMeasure{Aμ}, x::Ax) where {Aμ <: MappedArray, Ax <: AbstractArray}
#     μ.data
# end
