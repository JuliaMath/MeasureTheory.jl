export NormalInvChiSq
import MeasureBase
using MeasureTheory
using KeywordCalls

using SpecialFunctions: loggamma

# @parameterized NormalInvChiSq(μ, σ², κ, ν)

struct NormalInvChiSq{N,T} <: ParameterizedMeasure{N}
    par::NamedTuple{N,T}
    function (NormalInvChiSq{N,T}(nt::NamedTuple{N,T}) where {N,T<:Tuple})
        @assert nt.σ² > 0
        @assert nt.κ > 0
        @assert nt.ν > 0
        new{N,T}(nt)
    end
end

function NormalInvChiSq(μ, σ², κ, ν)
    NormalInvChiSq(NamedTuple{(:μ, :σ², :κ, :ν)}((μ, σ², κ, ν)))
end

@kwstruct NormalInvChiSq(μ, σ², κ, ν)

function Base.rand(
    ::ConstantRNG,
    ::Type{T},
    ::NormalInvChiSq{(:μ, :σ², :κ, :ν)},
) where {T<:Real}
    (; μ = zero(T), σ = one(T))
end

function MeasureBase.insupport(
    d::NormalInvChiSq,
    x::NamedTuple{(:μ, :σ),Tuple{M,S}},
) where {M<:Real,S<:Real}
    x.σ > 0
end

# See Gelman et al BDA 2nd edition, section 3.3
function MeasureBase.logdensity_def(
    d::NormalInvChiSq{(:μ, :σ², :κ, :ν)},
    x::NamedTuple{(:μ, :σ)},
)
    (μ₀, σ²₀, κ₀, ν₀) = params(d)

    (μ, σ) = x
    σ² = σ^2

    # TODO: Check the (μ,σ²) parameterization of a Normal, maybe change to that to avoid the sqrt.
    logdensity_def(ScaledInvChiSq(ν₀, σ²₀), σ²) +
    logdensity_def(Normal((μ = μ₀, σ² = σ² / κ₀)), μ)
end

function MeasureBase.basemeasure(::NormalInvChiSq)
    productmeasure((μ = LebesgueBase(), σ = LebesgueBase()))
end

# See Gelman et al BDA 2nd edition, section 3.3
function Base.rand(
    rng::AbstractRNG,
    ::Type{T},
    d::NormalInvChiSq{(:μ, :σ², :κ, :ν)},
) where {T}
    (μ₀, σ²₀, κ₀, ν₀) = params(d)

    σ² = rand(rng, T, ScaledInvChiSq(ν₀, σ²₀))
    σ = sqrt(σ²)
    # TODO: See about avoiding the sqrt here.
    μ = rand(rng, T, Normal(μ₀, σ / sqrt(κ₀)))
    return (; μ, σ)
end

function MeasureBase.pushfwd(::Type{Normal}, dist::NormalInvChiSq{(:μ, :σ², :κ, :ν)})
    (μ, σ², κ, ν) = params(dist)
    σ = sqrt(σ² * (1 + inv(κ)))
    return StudentT(ν, μ, σ)
end
