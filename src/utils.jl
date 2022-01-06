using LinearAlgebra

function solve(A::Union{AbstractMatrix, Factorization}, y::AbstractArray)
    (m, n) = size(A)
    n == 1 && return [dot(A, y) / sum(a -> a^2, A)]
    return A \ y
end
