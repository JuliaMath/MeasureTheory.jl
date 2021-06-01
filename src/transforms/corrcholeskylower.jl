
####
#### correlation cholesky Lower
####
# Adapted from https://github.com/tpapp/TransformVariables.jl/blob/master/src/special_arrays.jl

const TV = TransformVariables
using LinearAlgebra

"""
    CorrCholeskyLower(n)
Cholesky Lower of a correlation matrix of size `n`.
Transforms ``n×(n-1)/2`` real numbers to an ``n×n`` lower-triangular matrix `L`, such that
`L'*L` is a correlation matrix (positive definite, with unit diagonal).
# Notes
If
- `z` is a vector of `n` IID standard normal variates,
- `σ` is an `n`-element vector of standard deviations,
- `L` is obtained from `CorrCholeskyLower(n)`,
then `Diagonal(σ) * L' * z` will be a multivariate normal with the given variances and
correlation matrix `L' * L`.
"""
struct CorrCholeskyLower <: TV.VectorTransform
    n::Int
    function CorrCholeskyLower(n)
        @argcheck n ≥ 1 "Dimension should be positive."
        new(n)
    end
end

TV.dimension(t::CorrCholeskyLower) = TV.unit_triangular_dimension(t.n)

function TV.transform_with(flag::LogJacFlag, t::CorrCholeskyLower, x::AbstractVector, index)
    @unpack n = t
    T = extended_eltype(x)
    ℓ = logjac_zero(flag, T)
    L = Matrix{T}(undef, n, n)
    @inbounds for col in 1:n
        r = one(T)
        for row in col:n
            xi = x[index]
            L[row, col], r, ℓi = TV.l2_remainder_transform(flag, xi, r)
            ℓ += ℓi
            index += 1
        end
        L[col, col] = √r
    end
    # We could change the above code, but `transpose` is very cheap
    LowerTriangular(L), ℓ, index
end

TV.inverse_eltype(t::CorrCholeskyLower, L::LowerTriangular) = extended_eltype(L)

function TV.inverse_at!(x::AbstractVector, index, t::CorrCholeskyLower, L::LowerTriangular)
    @unpack n = t
    @argcheck size(L, 1) == n
    @inbounds for col in 1:n
        r = one(eltype(L))
        for row in col:n
            x[index], r = l2_remainder_inverse(L[row, col], r)
            index += 1
        end
    end
    index
end
