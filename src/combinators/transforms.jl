using TransformVariables: AbstractTransform, CallableTransform, CallableInverse

export Pushforward
export Pullback

struct Pushforward{F,M,L} <: AbstractMeasure
    f::F
    μ::M
    logjac::L
end

insupport(d::Pushforward, x) = insupport(d.μ, inverse(d.f)(x))

Pushforward(f, μ) = Pushforward(f, μ, True())

function Pretty.tile(pf::Pushforward{<:TV.CallableTransform})
    Pretty.list_layout(Pretty.tile.([pf.f.x, pf.μ, pf.logjac]); prefix = :Pushforward)
end

function Pretty.tile(pf::Pushforward)
    Pretty.list_layout(Pretty.tile.([pf.f, pf.μ, pf.logjac]); prefix = :Pushforward)
end
struct Pullback{F,M,L} <: AbstractMeasure
    f::F
    ν::M
    logjac::L
end

Pullback(f, ν) = Pullback(f, ν, True())

insupport(d::Pullback, x) = insupport(d.ν, d.f(x))

function Pretty.tile(pf::Pullback{<:TV.CallableTransform})
    Pretty.list_layout(Pretty.tile.([pf.f.x, pf.ν, pf.logjac]); prefix = :Pullback)
end

function Pretty.tile(pf::Pullback)
    Pretty.list_layout(Pretty.tile.([pf.f, pf.ν, pf.logjac]); prefix = :Pullback)
end

@inline function logdensity_def(pb::Pullback{F,M,True}, x) where {F<:CallableTransform,M}
    f = pb.f
    ν = pb.ν
    y, logJ = TV.transform_and_logjac(f.x, x)
    return logdensity_def(ν, y) + logJ
end

@inline function logdensity_def(pb::Pullback{F,M,False}, x) where {F<:CallableTransform,M}
    f = pb.f
    ν = pb.ν
    y = f(x)
    return logdensity_def(ν, y)
end

@inline function logdensity_def(pf::Pushforward{F,M,True}, y) where {F<:CallableTransform,M}
    f = pf.f
    μ = pf.μ
    x = TV.inverse(f.x)(y)
    _, logJ = TV.transform_and_logjac(f.x, x)
    return logdensity_def(μ, x) - logJ
end

@inline function logdensity_def(
    pf::Pushforward{F,M,False},
    y,
) where {F<:CallableTransform,M}
    f = pf.f
    μ = pf.μ
    x = TV.inverse(f.x)(y)
    return logdensity_def(μ, x)
end

Pullback(f::AbstractTransform, ν, logjac = True()) = Pullback(TV.transform(f), ν, logjac)

function Pushforward(f::AbstractTransform, ν, logjac = True())
    Pushforward(TV.transform(f), ν, logjac)
end

Pullback(f::CallableInverse, ν, logjac = True()) = Pushforward(TV.transform(f.x), ν, logjac)

Pushforward(f::CallableInverse, ν, logjac = True()) = Pullback(TV.transform(f.x), ν, logjac)

Base.rand(rng::AbstractRNG, T::Type, ν::Pushforward) = ν.f(rand(rng, T, ν.μ))

Base.rand(rng::AbstractRNG, T::Type, μ::Pullback) = μ.f(rand(rng, T, μ.ν))

testvalue(ν::Pushforward) = TV.transform(ν.f, testvalue(ν.μ))

testvalue(μ::Pullback) = TV.transform(TV.inverse(μ.f), testvalue(μ.ν))

basemeasure(μ::Pullback) = Pullback(μ.f, basemeasure(μ.ν), False())

basemeasure(ν::Pushforward) = Pushforward(ν.f, basemeasure(ν.μ), False())

as(ν::Pushforward) = ν.f ∘ as(ν.μ)

as(μ::Pullback) = TV.inverse(μ.f) ∘ μ.ν

as(::Lebesgue) = asℝ

# TODO: Make this work for affine embeddings
as(d::AffinePushfwd) = _as_affine(_firstval(d))

_firstval(d::AffinePushfwd) = first(values(getfield(getfield(d, :f), :par)))
_as_affine(x::Real) = asℝ
_as_affine(x::AbstractArray) = as(Vector, size(x, 1))

function basemeasure(
    ::Pushforward{TV.CallableTransform{T},Lebesgue{ℝ}},
) where {T<:TV.ScalarTransform}
    Lebesgue(ℝ)
end
function basemeasure(
    ::Pullback{TV.CallableTransform{T},Lebesgue{ℝ}},
) where {T<:TV.ScalarTransform}
    Lebesgue(ℝ)
end
# t = as𝕀
# μ = Normal()
# ν = Pushforward(t, μ)
# x = rand(μ)
# julia> logdensity_def(μ, x) ≈ logdensity_def(Pushforward(inverse(t), ν), x)
# true
