# Modified from
# https://github.com/tpapp/AltDistributions.jl

export LKJCholesky
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
struct LKJCholesky{k, N, T} <: ParameterizedMeasure{N}
    par :: NamedTuple{N,T}
end

LKJCholesky{k}(;kw...) where {k} = LKJCholesky{k}(NamedTuple(kw))

LKJCholesky{k}(nt::NamedTuple{N,T}) where {k,N,T} = LKJCholesky{k,N,T}(nt)

LKJCholesky(k,η) = LKJCholesky{k}(η)

LKJCholesky(k::Integer) = LKJCholesky(k, 1)

LKJCholesky{k}(η) where {k} = LKJCholesky{k, (:η,), Tuple{typeof(η)}}((η=η,))

function Base.show(io::IO, d::LKJCholesky{k}) where {k}
    println(io, "LKJCholesky{$k}(", getfield(d, :par), ")")
end

asparams(::Type{<:LKJCholesky}, ::Val{:η}) = asℝ₊
asparams(::Type{<:LKJCholesky}, ::Val{:logη}) = asℝ

# Modified from
# https://github.com/tpapp/AltDistributions.jl


using LinearAlgebra
using Tullio


function logdensity(d::LKJCholesky{k, (:η,), T}, L::Union{LinearAlgebra.AbstractTriangular, Diagonal}) where {k,T}
    η = d.η
    # z = diag(L)
    # sum(log.(z) .* ((k:-1:1) .+ 2*(η-1)))

    # Note: https://github.com/cscherrer/MeasureTheory.jl/issues/100#issuecomment-852428192
    c = k + 2(η - 1)
    @tullio s = (c - i) * log(L[i,i])
    return s
end

function logdensity(d::LKJCholesky{k, (:logη,), T}, L::Union{LinearAlgebra.AbstractTriangular, Diagonal}) where {k,T}
    c = k + 2 * expm1(d.logη)
    @tullio s = (c - i) * log(L[i,i])
    return s
end

asparams(::Type{<:LKJCholesky}, ::Val{:logη}) = asℝ₊

using TransformVariables

TransformVariables.as(::LKJCholesky{k}) where {k} = CorrCholesky(k)

# Should satisfy
# logdensity(basemeasure(d), rand(d)) == 0.0
function basemeasure(μ::LKJCholesky{k}) where {k}
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
function Base.rand(rng::AbstractRNG, ::Type, d::LKJCholesky{k, (:η,)}) where {k}
    return cholesky(Positive, rand(rng, Dists.LKJ(k, d.η)))
end;

function Base.rand(rng::AbstractRNG, ::Type, d::LKJCholesky{k, (:logη,)}) where {k}
    return cholesky(Positive, rand(rng, Dists.LKJ(k, exp(d.logη))))
end;

constructor(::Type{L}) where {k,L<:LKJCholesky{k}} = LKJCholesky{k}
