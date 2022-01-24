using LinearAlgebra
export mydot
using LazyArrays

function solve(A::Union{AbstractMatrix, Factorization}, y::AbstractArray)
    (m, n) = size(A)
    n == 1 && return [dot(A, y) / sum(a -> a^2, A)]
    return A \ y
end


@inline function mydot(a::Real, b::Real)
    return a * b
end

@inline function mydot(a::Ta, b::Tb) where {Ta<:AbstractArray, Tb<:AbstractArray}
    return dot(a, b)
end

@inline function _mydot(a::NTuple{0}, b::NTuple{0}, s) 
    return s
end

@inline function mydot(a::Tuple, b::Tuple) 
    z = 0.0
    _mydot(a, b, z)
end

@inline function _mydot(a::Tuple, b::Tuple, s)
    _mydot(Base.tail(a), Base.tail(b), s + mydot(first(a), first(b)))
end

@generated function mydot(a::SVector{N,Ta}, b::SVector{N,Tb}) where {N,Ta,Tb}
    z = zero(float(Base.promote_eltype(Ta, Tb)))
    quote
        $(Expr(:meta, :inline))
        result = $z
        @inbounds Base.Cartesian.@nexprs $N i -> begin
            result += a[i] * b[i]
        end
        return result
    end
end

in𝕀(x) = static(0.0) ≤ x ≤ static(1.0)
inℝ₊(x) = static(0.0) ≤ x