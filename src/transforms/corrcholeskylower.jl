# Adapted from https://github.com/tpapp/TransformVariables.jl/blob/master/src/special_arrays.jl

using LinearAlgebra
const CorrCholeskyUpper = CorrCholeskyFactor

export CorrCholeskyLower

"""
    CorrCholeskyLower(n)
Lower Cholesky factor of a correlation matrix of size `n`.
Transforms ``n(n-1)/2`` real numbers to an ``n×n`` lower-triangular matrix `L`, such that
`L*L'` is a correlation matrix (positive definite, with unit diagonal).

# Notes
If
- `z` is a vector of `n` IID standard normal variates,
- `σ` is an `n`-element vector of standard deviations,
- `L` is obtained from `CorrCholeskyLower(n)`,
then `Diagonal(σ) * L * z` is a zero-centered multivariate normal variate with the standard deviations `σ` and
correlation matrix `L * L'`.
"""
struct CorrCholeskyLower <: TV.VectorTransform
    n::Int
end

TV.dimension(t::CorrCholeskyLower) = TV.dimension(CorrCholeskyUpper(t.n))

# TODO: For now we just transpose the CorrCholeskyUpper result. We should
# consider whether it can help performance to implement this directly for the
# lower triangular case
function TV.transform_with(flag::TV.LogJacFlag, t::CorrCholeskyLower, x::AbstractVector, index)
    U, ℓ, index = TV.transform_with(flag, CorrCholeskyUpper(t.n), x, index)
    return U', ℓ, index
end

TV.inverse_eltype(t::CorrCholeskyLower, L::LowerTriangular) = TV.extended_eltype(L)

function TV.inverse_at!(x::AbstractVector, index, t::CorrCholeskyLower, L::LowerTriangular)
    return TV.inverse_at!(x, index, CorrCholeskyUpper(t.n), L')
end
