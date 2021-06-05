# Modified from
# https://github.com/tpapp/AltDistributions.jl

export LKJL
using PositiveFactorizations

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

LKJL{k}(;kw...) where {k} = LKJL{k}(NamedTuple(kw))

LKJL{k}(nt::NamedTuple{N,T}) where {k,N,T} = LKJL{k,N,T}(nt)

LKJL(k,η) = LKJL{k}(η)

LKJL(k::Integer) = LKJL(k, 1)

LKJL{k}(η) where {k} = LKJL{k, (:η,), Tuple{typeof(η)}}((η=η,))

function Base.show(io::IO, d::LKJL{k}) where {k}
    println(io, "LKJL{$k}(", getfield(d, :par), ")")
end

asparams(::Type{<:LKJL}, ::Val{:η}) = asℝ₊
asparams(::Type{<:LKJL}, ::Val{:logη}) = asℝ

# Modified from
# https://github.com/tpapp/AltDistributions.jl


using LinearAlgebra
using Tullio


function logdensity(d::LKJL{k, (:η,), T}, L::Union{LinearAlgebra.AbstractTriangular, Diagonal}) where {k,T}
    η = d.η
    # z = diag(L)
    # sum(log.(z) .* ((k:-1:1) .+ 2*(η-1)))

    # Note: https://github.com/cscherrer/MeasureTheory.jl/issues/100#issuecomment-852428192
    c = k + 2(η - 1)
    @tullio s = (c - i) * log(L[i,i])
    return s
end

function logdensity(d::LKJL{k, (:logη,), T}, L::Union{LinearAlgebra.AbstractTriangular, Diagonal}) where {k,T}
    c = k + 2 * expm1(d.logη)
    @tullio s = (c - i) * log(L[i,i])
    return s
end

asparams(::Type{<:LKJL}, ::Val{:logη}) = asℝ₊

using TransformVariables

TransformVariables.as(::LKJL{k}) where {k} = CorrCholeskyLower(k)

# Should satisfy
# logdensity(basemeasure(d), rand(d)) == 0.0
function basemeasure(μ::LKJL{k}) where {k}
    t = as(μ)
    d = dimension(t)
    return Pushforward(t, Lebesgue(ℝ)^d, false)
end

# Note (from @sethaxen): this can be done without the cholesky decomposition, by randomly
# sampling the canonical partial correlations, as in the LKJ paper, and then
# mapping them to the cholesky factor instead of the correlation matrix. Stan
# implements* this. But that can be for a future PR. 
#
# * https://github.com/stan-dev/math/blob/develop/stan/math/prim/prob/lkj_corr_cholesky_rng.hpp
function Base.rand(rng::AbstractRNG, ::Type, d::LKJL{k, (:η,)}) where {k}
    return cholesky(rand(rng, Dists.LKJ(k, d.η))).L
end;

function Base.rand(rng::AbstractRNG, ::Type, d::LKJL{k, (:logη,)}) where {k}
    return cholesky(Positive, rand(rng, Dists.LKJ(k, exp(d.logη)))).L
end;

constructor(::Type{L}) where {k,L<:LKJL{k}} = LKJL{k}
