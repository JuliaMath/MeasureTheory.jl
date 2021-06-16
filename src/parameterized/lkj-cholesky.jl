# Modified from
# https://github.com/tpapp/AltDistributions.jl

# See also the Stan manual (the "Stanual", though nobody calls it that)
# https://mc-stan.org/docs/2_27/reference-manual/cholesky-factors-of-correlation-matrices-1.html

export LKJCholesky
using PositiveFactorizations

const TV = TransformVariables

const CorrCholeskyUpper = TV.CorrCholeskyFactor

"""
    LKJCholesky{k}(η=1.0)
    LKJCholesky{k}(logη=0.0)

`LKJCholesky{k}` gives the `k×k` LKJ distribution (Lewandowski et al 2009)
expressed as a Cholesky decomposition. As a special case, for
`C = rand(LKJCholesky{k}(η=1.0))` (or equivalently
`C=rand(LKJCholesky{k}(logη=0.0))`), `C.L * C.U` is uniform over the set of all
`k×k` correlation matrices. NOte, however, that in this case `C.L` and `C.U` are
**not** sampled uniformly (because the multiplication is nonlinear).

The `logdensity` method for this measure applies for `LowerTriangular`,
`UpperTriangular`, or `Diagonal` matrices, and will "do the right thing". The
`logdensity` **does not check if `L*U` yields a valid correlation matrix**.

Valid values are ``η > 0``. When ``η > 1``, the distribution is unimodal with a
peak at `I`, while ``0 < η < 1`` yields a trough. ``η = 2`` is recommended as a vague prior.

Adapted from
https://github.com/tpapp/AltDistributions.jl
"""
struct LKJCholesky{k, N, T} <: ParameterizedMeasure{N}
    par :: NamedTuple{N,T}

    function LKJCholesky{k,N,T}(nt) where {k,N,T}
        @warn """
        WORK IN PROGRESS
        
        `LKJCholesky` does not yet have the correct base measure
        """

        new{k,N,T}(nt)
    end
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

logdensity(d::LKJCholesky, C::Cholesky) = logdensity(d, C.UL)

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


TV.as(::LKJCholesky{k}) where {k} = CorrCholesky(k)

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
