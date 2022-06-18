export Affine, AffineTransform
using LinearAlgebra
import Base

const AFFINEPARS = [
    (:μ, :σ)
    (:μ, :λ)
    (:σ,)
    (:λ,)
    (:μ,)
]

struct AffineTransform{N,T}
    par::NamedTuple{N,T}
end

function Pretty.tile(f::AffineTransform)
    result = Pretty.literal("AffineTransform")
    result *= Pretty.literal(sprint(show, params(f); context = :compact => true))
    result
end

Base.show(io::IO, f::AffineTransform) = Pretty.pprint(io, f)

params(f::AffineTransform) = getfield(f, :par)

@inline Base.getproperty(d::AffineTransform, s::Symbol) = getfield(getfield(d, :par), s)

Base.propertynames(d::AffineTransform{N}) where {N} = N

import InverseFunctions: inverse

@inline inverse(f::AffineTransform{(:μ, :σ)}) = AffineTransform((μ = -(f.σ \ f.μ), λ = f.σ))
@inline inverse(f::AffineTransform{(:μ, :λ)}) = AffineTransform((μ = -f.λ * f.μ, σ = f.λ))
@inline inverse(f::AffineTransform{(:σ,)}) = AffineTransform((λ = f.σ,))
@inline inverse(f::AffineTransform{(:λ,)}) = AffineTransform((σ = f.λ,))
@inline inverse(f::AffineTransform{(:μ,)}) = AffineTransform((μ = -f.μ,))

# `size(f) == (m,n)` means `f : ℝⁿ → ℝᵐ`  
Base.size(f::AffineTransform{(:μ, :σ)}) = size(f.σ)
Base.size(f::AffineTransform{(:μ, :λ)}) = size(f.λ)
Base.size(f::AffineTransform{(:σ,)}) = size(f.σ)
Base.size(f::AffineTransform{(:λ,)}) = size(f.λ)

function Base.size(f::AffineTransform{(:μ,)})
    (n,) = size(f.μ)
    return (n, n)
end

Base.size(f::AffineTransform, n::Int) = @inbounds size(f)[n]

(f::AffineTransform{(:μ,)})(x) = x + f.μ
(f::AffineTransform{(:σ,)})(x) = f.σ * x
(f::AffineTransform{(:λ,)})(x) = f.λ \ x
(f::AffineTransform{(:μ, :σ)})(x) = f.σ * x + f.μ
(f::AffineTransform{(:μ, :λ)})(x) = f.λ \ x + f.μ

rowsize(x) = ()
rowsize(x::AbstractArray) = (size(x, 1),)

function rowsize(f::AffineTransform)
    size_f = size(f)

    size_f isa Tuple{} && return 0
    return first(size_f)
end

colsize(x) = ()
colsize(x::AbstractArray) = (size(x, 2),)

function colsize(f::AffineTransform)
    size_f = size(f)

    size_f isa NTuple{2} && return last(size_f)
    return 0
end

@inline function apply!(x, f::AffineTransform{(:μ,)}, z)
    x .= z .+ f.μ
    return x
end

@inline function apply!(x, f::AffineTransform{(:σ,)}, z)
    mul!(x, f.σ, z)
    return x
end

@inline function apply!(x, f::AffineTransform{(:λ,),Tuple{F}}, z) where {F<:Factorization}
    ldiv!(x, f.λ, z)
    return x
end

@inline function apply!(x, f::AffineTransform{(:λ,)}, z)
    ldiv!(x, factorize(f.λ), z)
    return x
end

@inline function apply!(x, f::AffineTransform{(:μ, :σ)}, z)
    apply!(x, AffineTransform((σ = f.σ,)), z)
    apply!(x, AffineTransform((μ = f.μ,)), x)
    return x
end

@inline function apply!(x, f::AffineTransform{(:μ, :λ)}, z)
    apply!(x, AffineTransform((λ = f.λ,)), z)
    apply!(x, AffineTransform((μ = f.μ,)), x)
    return x
end

@inline function logjac(x::AbstractMatrix)
    (m, n) = size(x)
    m == n && return first(logabsdet(x))

    # Equivalent to sum(log, svdvals(x)), but much faster
    m > n && return first(logabsdet(x' * x)) / 2
    return first(logabsdet(x * x')) / 2
end

logjac(x::Number) = log(abs(x))

# TODO: `log` doesn't work for the multivariate case, we need the log absolute determinant
logjac(f::AffineTransform{(:μ, :σ)}) = logjac(f.σ)
logjac(f::AffineTransform{(:μ, :λ)}) = -logjac(f.λ)
logjac(f::AffineTransform{(:σ,)}) = logjac(f.σ)
logjac(f::AffineTransform{(:λ,)}) = -logjac(f.λ)
logjac(f::AffineTransform{(:μ,)}) = 0.0

###############################################################################

struct OrthoLebesgue{N,T} <: PrimitiveMeasure
    par::NamedTuple{N,T}

    OrthoLebesgue(nt::NamedTuple{N,T}) where {N,T} = new{N,T}(nt)
end

basemeasure(d::OrthoLebesgue) = d

logdensity_def(::OrthoLebesgue, x) = static(0)

struct Affine{N,M,T} <: AbstractMeasure
    f::AffineTransform{N,T}
    parent::M
end

function Pretty.tile(d::Affine)
    pars = Pretty.literal(sprint(show, params(d.f); context = :compact => true))

    Pretty.list_layout([pars, Pretty.tile(d.parent)]; prefix = :Affine)
end

function testvalue(d::Affine)
    f = getfield(d, :f)
    z = testvalue(parent(d))
    return f(z)
end

Affine(nt::NamedTuple, μ::AbstractMeasure) = affine(nt, μ)

Affine(nt::NamedTuple) = affine(nt)

Base.parent(d::Affine) = getfield(d, :parent)

function params(μ::Affine)
    nt1 = getfield(getfield(μ, :f), :par)
    nt2 = params(parent(μ))
    return merge(nt1, nt2)
end

function paramnames(::Type{A}) where {N,M,A<:Affine{N,M}}
    tuple(union(N, paramnames(M))...)
end

Base.propertynames(d::Affine{N}) where {N} = N ∪ (:parent, :f)

@inline function Base.getproperty(d::Affine, s::Symbol)
    if s === :parent
        return getfield(d, :parent)
    elseif s === :f
        return getfield(d, :f)
    else
        return getproperty(getfield(d, :f), s)
    end
end

Base.size(d::Affine) = size(d.μ)
Base.size(d::Affine{(:σ,)}) = (size(d.σ, 1),)
Base.size(d::Affine{(:λ,)}) = (size(d.λ, 2),)

@inline function logdensity_def(d::Affine{(:σ,)}, x::AbstractArray)
    z = solve(d.σ, x)
    MeasureBase.unsafe_logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:λ,)}, x::AbstractArray)
    z = d.λ * x
    MeasureBase.unsafe_logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:μ,)}, x::AbstractArray)
    z = mappedarray(-, x, d.μ)
    MeasureBase.unsafe_logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:μ, :σ)}, x::AbstractArray)
    z = d.σ \ mappedarray(-, x, d.μ)
    MeasureBase.unsafe_logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:μ, :λ)}, x::AbstractArray)
    z = d.λ * mappedarray(-, x, d.μ)
    MeasureBase.unsafe_logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine, x)
    z = inverse(d.f)(x)
    logdensity_def(d.parent, z)
end

# # # logdensity_def(d::Affine{(:μ,:λ)}, x) = logdensity_def(d.parent, d.σ \ (x - d.μ))
# # @inline function logdensity_def(d::Affine{(:μ,:σ), P, Tuple{V,M}}, x) where {P, V<:AbstractVector, M<:AbstractMatrix}
# #     z = x - d.μ
# #     σ = d.σ
# #     if σ isa Factorization
# #         ldiv!(σ, z)
# #     else
# #         ldiv!(factorize(σ), z)
# #     end
# #     sum(zⱼ -> logdensity_def(d.parent, zⱼ), z)
# # end

# # # logdensity_def(d::Affine{(:μ,:λ)}, x) = logdensity_def(d.parent, d.λ * (x - d.μ))
# # @inline function logdensity_def(d::Affine{(:μ,:λ), P,Tuple{V,M}}, x) where {P,V<:AbstractVector, M<:AbstractMatrix}
# #     z = x - d.μ
# #     lmul!(d.λ, z)
# #     logdensity_def(d.parent, z)
# # end

@inline function basemeasure(d::Affine{N,M,Tuple{A}}) where {N,M,A<:AbstractArray}
    weightedmeasure(-logjac(d), OrthoLebesgue(params(d)))
end

@inline function basemeasure(
    d::Affine{N,M,Tuple{A1,A2}},
) where {N,M,A1<:AbstractArray,A2<:AbstractArray}
    weightedmeasure(-logjac(d), OrthoLebesgue(params(d)))
end

@inline basemeasure(d::Affine) = affine(getfield(d, :f), basemeasure(d.parent))

# We can't do this until we know we're working with Lebesgue measure, since for
# example it wouldn't make sense to apply a log-Jacobian to a point measure
@inline function basemeasure(d::Affine{N,L}) where {N,L<:Lebesgue}
    weightedmeasure(-logjac(d), d.parent)
end
@inline function basemeasure(d::Affine{N,L}) where {N,L<:LebesgueMeasure}
    weightedmeasure(-logjac(d), d.parent)
end

logjac(d::Affine) = logjac(getfield(d, :f))

function Random.rand!(
    rng::Random.AbstractRNG,
    d::Affine,
    x::AbstractVector{T},
    z = Vector{T}(undef, size(getfield(d, :f), 2)),
) where {T}
    rand!(rng, parent(d), z)
    f = getfield(d, :f)
    apply!(x, f, z)
    return x
end

# function Base.rand(rng::Random.AbstractRNG, ::Type{T}, d::Affine) where {T}
#     f = getfield(d, :f)
#     z = rand(rng, T, parent(d))
#     apply!(x, f, z)
#     return z
# end

supportdim(nt::NamedTuple{(:μ, :σ)}) = colsize(nt.σ)
supportdim(nt::NamedTuple{(:μ, :λ)}) = rowsize(nt.λ)
supportdim(nt::NamedTuple{(:σ,)}) = colsize(nt.σ)
supportdim(nt::NamedTuple{(:λ,)}) = rowsize(nt.λ)
supportdim(nt::NamedTuple{(:μ,)}) = size(nt.μ)

asparams(::Affine, ::StaticSymbol{:μ}) = asℝ
asparams(::Affine, ::StaticSymbol{:σ}) = asℝ₊
asparams(::Type{A}, ::StaticSymbol{:μ}) where {A<:Affine} = asℝ
asparams(::Type{A}, ::StaticSymbol{:σ}) where {A<:Affine} = asℝ₊

function asparams(d::Affine{N,M,T}, ::StaticSymbol{:μ}) where {N,M,T<:AbstractArray}
    as(Array, asℝ, size(d.μ))
end

function asparams(d::Affine{N,M,T}, ::StaticSymbol{:σ}) where {N,M,T<:AbstractArray}
    as(Array, asℝ, size(d.σ))
end

function Base.rand(rng::Random.AbstractRNG, ::Type{T}, d::Affine) where {T}
    z = rand(rng, T, parent(d))
    f = getfield(d, :f)
    return f(z)
end

@inline function insupport(d::Affine, x)
    insupport(d.parent, inverse(d.f)(x))
end

@inline function Distributions.cdf(d::Affine, x) 
    cdf(parent(d), inverse(d.f)(x))
end