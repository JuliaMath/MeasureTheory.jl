using TransformVariables
using TransformVariables: AbstractTransform, InverseTransform

export Pushforward
export Pullback

struct Pushforward{F,M} <: AbstractMeasure
    f::F
    μ::M
end

struct Pullback{F,M} <: AbstractMeasure
    f::F
    ν::M
end

function logdensity(pb::Pullback{F}, x) where {F <: AbstractTransform}
    f = pb.f
    ν = pb.ν
    y, ℓ = transform_and_logjac(f, x)
    logdensity(ν,y) + ℓ
end

function logdensity(pf::Pushforward{F}, y) where {F <: AbstractTransform}
    f = pf.f
    μ = pf.μ
    x = inverse(f, y)
    _, ℓ = transform_and_logjac(f, x)
    logdensity(μ, x) - ℓ
end

Pullback(f::InverseTransform, ν) = Pushforward(f.transform, ν)

Pushforward(f::InverseTransform, ν) = Pullback(f.transform, ν)

Base.rand(rng::AbstractRNG, T::Type, ν::Pushforward) = ν.f(rand(rng, ν.μ))

Base.rand(rng::AbstractRNG, T::Type, μ::Pullback) = μ.f(rand(rng, μ.ν))

testvalue(ν::Pushforward) = transform(ν.f, testvalue(ν.μ))

testvalue(μ::Pullback) = transform(μ.f, testvalue(μ.ν))


# t = as𝕀
# μ = Normal()
# ν = Pushforward(t, μ)
# x = rand(μ)
# julia> logdensity(μ, x) ≈ logdensity(Pushforward(inverse(t), ν), x)
# true
