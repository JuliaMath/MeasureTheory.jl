export Affine, AffineTransform
using LinearAlgebra

const AFFINEPARS = [
    (:μ, :σ)
    (:μ, :ω)
    (:σ,)
    (:ω,)
    (:μ,)
]

struct AffineTransform{N,T}
    par::NamedTuple{N,T}
end

quoteof(f::AffineTransform) = :(AffineTransform($(quoteof(f.par))))

params(f::AffineTransform) = getfield(f, :par)

@inline Base.getproperty(d::AffineTransform, s::Symbol) = getfield(getfield(d, :par), s)

Base.propertynames(d::AffineTransform{N}) where {N} = N

import InverseFunctions: inverse

@inline inverse(f::AffineTransform{(:μ, :σ)}) =
    AffineTransform((μ = -(f.σ \ f.μ), ω = f.σ))
@inline inverse(f::AffineTransform{(:μ, :ω)}) = AffineTransform((μ = -f.ω * f.μ, σ = f.ω))
@inline inverse(f::AffineTransform{(:σ,)}) = AffineTransform((ω = f.σ,))
@inline inverse(f::AffineTransform{(:ω,)}) = AffineTransform((σ = f.ω,))
@inline inverse(f::AffineTransform{(:μ,)}) = AffineTransform((μ = -f.μ,))

# `size(f) == (m,n)` means `f : ℝⁿ → ℝᵐ`  
Base.size(f::AffineTransform{(:μ, :σ)}) = size(f.σ)
Base.size(f::AffineTransform{(:μ, :ω)}) = size(f.ω)
Base.size(f::AffineTransform{(:σ,)}) = size(f.σ)
Base.size(f::AffineTransform{(:ω,)}) = size(f.ω)

function Base.size(f::AffineTransform{(:μ,)})
    (n,) = size(f.μ)
    return (n, n)
end

Base.size(f::AffineTransform, n::Int) = @inbounds size(f)[n]

(f::AffineTransform{(:μ,)})(x) = x + f.μ
(f::AffineTransform{(:σ,)})(x) = f.σ * x
(f::AffineTransform{(:ω,)})(x) = f.ω \ x
(f::AffineTransform{(:μ, :σ)})(x) = f.σ * x + f.μ
(f::AffineTransform{(:μ, :ω)})(x) = f.ω \ x + f.μ

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

@inline function apply!(x, f::AffineTransform{(:ω,),Tuple{F}}, z) where {F<:Factorization}
    ldiv!(x, f.ω, z)
    return x
end

@inline function apply!(x, f::AffineTransform{(:ω,)}, z)
    ldiv!(x, factorize(f.ω), z)
    return x
end

@inline function apply!(x, f::AffineTransform{(:μ, :σ)}, z)
    apply!(x, AffineTransform((σ = f.σ,)), z)
    apply!(x, AffineTransform((μ = f.μ,)), x)
    return x
end

@inline function apply!(x, f::AffineTransform{(:μ, :ω)}, z)
    apply!(x, AffineTransform((ω = f.ω,)), z)
    apply!(x, AffineTransform((μ = f.μ,)), x)
    return x
end

function logjac(x::AbstractMatrix)
    (m, n) = size(x)
    m == n && return first(logabsdet(x))

    # Equivalent to sum(log, svdvals(x)), but much faster
    m > n && return first(logabsdet(x' * x)) / 2
    return first(logabsdet(x * x')) / 2
end

logjac(x::Number) = log(abs(x))

# TODO: `log` doesn't work for the multivariate case, we need the log absolute determinant
logjac(f::AffineTransform{(:μ, :σ)}) = logjac(f.σ)
logjac(f::AffineTransform{(:μ, :ω)}) = -logjac(f.ω)
logjac(f::AffineTransform{(:σ,)}) = logjac(f.σ)
logjac(f::AffineTransform{(:ω,)}) = -logjac(f.ω)
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
    Pretty.list_layout(Pretty.tile.([params(d.f), d.parent]); prefix=:Affine)
end

function testvalue(d::Affine)
    f = getfield(d, :f)
    z = testvalue(parent(d))
    return f(z)
end

Affine(nt::NamedTuple, μ::AbstractMeasure) = affine(nt, μ)

Affine(nt::NamedTuple) = affine(nt)

parent(d::Affine) = getfield(d, :parent)

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
Base.size(d::Affine{(:ω,)}) = (size(d.ω, 2),)

@inline function logdensity_def(d::Affine{(:σ,)}, x::AbstractArray)
    z = solve(d.σ, x)
    MeasureBase._logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:ω,)}, x::AbstractArray)
    z = d.ω * x
    MeasureBase._logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:μ,)}, x::AbstractArray)
    z = mappedarray(-, x, d.μ)
    MeasureBase._logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:μ, :σ)}, x::AbstractArray)
    z = d.σ \ mappedarray(-, x, d.μ)
    MeasureBase._logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine{(:μ, :ω)}, x::AbstractArray)
    z = d.ω * mappedarray(-, x, d.μ)
    MeasureBase._logdensityof(d.parent, z)
end

@inline function logdensity_def(d::Affine, x)
    z = inverse(d.f)(x)
    logdensity_def(d.parent, z)
end

# # # logdensity_def(d::Affine{(:μ,:ω)}, x) = logdensity_def(d.parent, d.σ \ (x - d.μ))
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

# # # logdensity_def(d::Affine{(:μ,:ω)}, x) = logdensity_def(d.parent, d.ω * (x - d.μ))
# # @inline function logdensity_def(d::Affine{(:μ,:ω), P,Tuple{V,M}}, x) where {P,V<:AbstractVector, M<:AbstractMatrix}
# #     z = x - d.μ
# #     lmul!(d.ω, z)
# #     logdensity_def(d.parent, z)
# # end
 
@inline basemeasure(d::Affine{N,M,Tuple{Vararg{T,K}}}) where {K,N,M,T<:AbstractArray} = affine(getfield(d, :f), rootmeasure(d.parent))

@inline basemeasure(d::Affine) = affine(getfield(d, :f), basemeasure(d.parent))

# We can't do this until we know we're working with Lebesgue measure, since for
# example it wouldn't make sense to apply a log-Jacobian to a point measure
@inline basemeasure(d::Affine{N,L}) where {N,L<:Lebesgue} = weightedmeasure(-logjac(d), d.parent)
@inline basemeasure(d::Affine{N,L}) where {N,L<:LebesgueMeasure} = weightedmeasure(-logjac(d), d.parent)


@inline function basemeasure(d::Affine{N,P, Tuple{Vararg{T, K}}}) where {K,N,L<:Union{<:Lebesgue, <:LebesgueMeasure},P<:PowerMeasure{L}, T<:AbstractArray} 
    weightedmeasure(-logjac(d), OrthoLebesgue(params(d)))
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
supportdim(nt::NamedTuple{(:μ, :ω)}) = rowsize(nt.ω)
supportdim(nt::NamedTuple{(:σ,)}) = colsize(nt.σ)
supportdim(nt::NamedTuple{(:ω,)}) = rowsize(nt.ω)
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
