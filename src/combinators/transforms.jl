using TransformVariables
using TransformVariables: AbstractTransform, CallableTransform, CallableInverse

export Pushforward
export Pullback

struct Pushforward{F,M} <: AbstractMeasure
    f::F
    μ::M
    logjac::Bool
end

Pushforward(f, μ) = Pushforward(f, μ, true)

struct Pullback{F,M} <: AbstractMeasure
    f::F
    ν::M
    logjac::Bool
end

Pullback(f, ν) = Pullback(f, ν, true)

@inline function logdensity_def(pb::Pullback{F}, x) where {F<:CallableTransform}
    f = pb.f
    ν = pb.ν
    if pb.logjac
        y, logJ = transform_and_logjac(f.t, x)
        return logdensity_def(ν, y) + logJ
    else
        y = f(x)
        return logdensity_def(ν, y)
    end
end

@inline function logdensity_def(pf::Pushforward{F}, y) where {F<:CallableTransform}
    f = pf.f
    μ = pf.μ
    x = inverse(f.t)(y)
    if pf.logjac
        _, logJ = transform_and_logjac(f.t, x)
        return logdensity_def(μ, x) - logJ
    else
        return logdensity_def(μ, x)
    end
end

Pullback(f::AbstractTransform, ν, logjac::Bool = true) = Pullback(transform(f), ν, logjac)
Pushforward(f::AbstractTransform, ν, logjac::Bool = true) =
    Pushforward(transform(f), ν, logjac)

Pullback(f::CallableInverse, ν, logjac::Bool = true) =
    Pushforward(transform(f.t), ν, logjac)

Pushforward(f::CallableInverse, ν, logjac::Bool = true) =
    Pullback(transform(f.t), ν, logjac)

Base.rand(rng::AbstractRNG, T::Type, ν::Pushforward) = ν.f(rand(rng, T, ν.μ))

Base.rand(rng::AbstractRNG, T::Type, μ::Pullback) = μ.f(rand(rng, T, μ.ν))

testvalue(ν::Pushforward) = transform(ν.f, testvalue(ν.μ))

testvalue(μ::Pullback) = transform(inverse(μ.f), testvalue(μ.ν))

basemeasure(μ::Pullback) = Pullback(μ.f, basemeasure(μ.ν), false)

basemeasure(ν::Pushforward) = Pushforward(ν.f, basemeasure(ν.μ), false)

TV.as(ν::Pushforward) = ν.f ∘ as(ν.μ)

TV.as(μ::Pullback) = inverse(μ.f) ∘ μ.ν

TV.as(::Lebesgue) = asℝ

# TODO: Make this work for affine embeddings
TV.as(d::Affine) = _as_affine(_firstval(d))

_firstval(d::Affine) = first(values(getfield(getfield(d, :f), :par)))
_as_affine(x::Real) = asℝ
_as_affine(x::AbstractArray) = as(Vector, size(x, 1))

basemeasure(
    ::Pushforward{TV.CallableTransform{T},Lebesgue{ℝ}},
) where {T<:TV.ScalarTransform} = Lebesgue(ℝ)
basemeasure(::Pullback{TV.CallableTransform{T},Lebesgue{ℝ}}) where {T<:TV.ScalarTransform} =
    Lebesgue(ℝ)
# t = as𝕀
# μ = Normal()
# ν = Pushforward(t, μ)
# x = rand(μ)
# julia> logdensity_def(μ, x) ≈ logdensity_def(Pushforward(inverse(t), ν), x)
# true
