using LinearAlgebra

function solve(A::Union{AbstractMatrix, Factorization}, y::AbstractArray)
    (m, n) = size(A)
    n == 1 && return [dot(A, y) / sum(a -> a^2, A)]
    return A \ y
end


@generated function mydot(a::SVector{N,Ta}, b::SVector{N,Tb}) where {N,Ta,Tb}
    z = zero(promote_type(Ta, Tb))
    quote
        $(Expr(:meta, :inline))
        result = $z
        @inbounds Base.Cartesian.@nexprs $N i -> begin
            result += a[i] * b[i]
        end
        return result
    end
end
