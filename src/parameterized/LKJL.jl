# Modified from
# https://github.com/tpapp/AltDistributions.jl

export LKJL

"""
The LKJ distribution (Lewandowski et al 2009) for the Cholesky factor L of correlation
matrices.

A correlation matrix ``Ω=L*L'`` has a density proportional to ``|Ω|^{η-1}``. However, it is usually not
necessary to construct ``Ω``, so this distribution is formulated for the Cholesky
decomposition `L*L'` and is supported on the lower triangular cholesky factor `L`.

Note that the method **does not check if `L` yields a valid correlation matrix**.
Valid values are ``η > 0``. When ``η > 1``, the distribution is unimodal at `Ω=I`, while
``0 < η < 1`` has a trough. ``η = 2`` is recommended as a vague prior.
When ``η = 1``, the density is uniform in `Ω`, but not in `L`, because of the Jacobian
correction of the transformation.

Adapted from
https://github.com/tpapp/AltDistributions.jl
"""
struct LKJL{k, N, T} <: ParameterizedMeasure{N}
    par :: NamedTuple{N,T}
end

LKJL{k}(nt::NamedTuple{N,T}) where {k,N,T} = LKJL{k,N,T}(nt)

LKJL(k,η) = LKJL{k}(η)

LKJL(k) = LKJL(k, 1.0)

LKJL{k}(η) where {k} = LKJL{k, (:η,), Tuple{typeof(η)}}((η=η,))

function Base.show(io::IO, d::LKJL{k}) where {k}
    println(io, "LKJL{$k}(", getfield(d, :par), ")")
end

asparams(::Type{<:LKJL}, ::Val{:η}) = asℝ₊

# Modified from
# https://github.com/tpapp/AltDistributions.jl


using LinearAlgebra
using Tullio


function logdensity(d::LKJL{k}, L::Union{LinearAlgebra.AbstractTriangular, Diagonal}) where {k}
    η = d.η
    # z = diag(L)
    # sum(log.(z) .* ((k:-1:1) .+ 2*(η-1)))

    # Note: https://github.com/cscherrer/MeasureTheory.jl/issues/100#issuecomment-852428192
    c = k + 2(η - 1)
    @tullio s = (c - i) * log(L[i,i])
    return s
end

using TransformVariables

TransformVariables.as(::LKJL{k}) where {k} = CorrCholeskyLower(k)

# TODO: I think this is wrong
function basemeasure(μ::LKJL{k}) where {k}
    t = as(μ)
    d = dimension(t)
    return Pushforward(t, Lebesgue(ℝ)^d, false)
end

function Base.rand(rng::AbstractRNG, ::Type, d::LKJL{k}) where {k}
    return cholesky(rand(rng, Dists.LKJ(k, d.η))).L
end;
