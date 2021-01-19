import Base
import FillArrays

PowerMeasure{M,N,D} = ProductMeasure{Fill{M,N,D}}

function Base.show(io::IO, μ::PowerMeasure)
    io = IOContext(io, :compact => true)
    print(io, μ.data.value, " ^ ", size(μ.data))
end

function Base.show_unquoted(io::IO, μ::PowerMeasure{M,N,D}, indent::Int, prec::Int) where {M,N,D}
    io = IOContext(io, :compact => true)
    if Base.operator_precedence(:^) ≤ prec
        print(io, "(")
        show(io, μ.data.value)
        print(io, ")")
    else
        show(io, size(μ.data))
    end
    return nothing
end

Base.:^(μ::AbstractMeasure, n::Integer) = μ ^ (n,)

Base.:^(μ::AbstractMeasure, size::NTuple{N,I}) where {N, I <: Integer} = ProductMeasure(Fill(μ, size))

sampletype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{sampletype(d.data[1]), N}

function Base.:^(μ::WeightedMeasure, n::NTuple{N,Int}) where {N}
    k = prod(n) * μ.logweight
    return WeightedMeasure(k, μ.base^n)
end

basemeasure(μ::PowerMeasure) = basemeasure(μ.data[1])^size(μ.data)
