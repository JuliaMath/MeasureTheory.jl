using TransformVariables
using TransformVariables: AbstractTransform, InverseTransform

export Pushforward
export Pullback

struct Pushforward{F,M} <: AbstractMeasure
    f::F
    μ::M
    logjac::Bool
end

Pushforward(f,μ) = Pushforward(f, μ, true)

struct Pullback{F,M} <: AbstractMeasure
    f::F
    ν::M
    logjac::Bool
end

Pullback(f,ν) = Pullback(f, ν, true)

function logdensity(pb::Pullback{F}, x) where {F <: AbstractTransform}
    f = pb.f
    ν = pb.ν
    if pb.logjac
        y, ℓ = transform_and_logjac(f, x)
        return logdensity(ν, y) + ℓ
    else
        y = transform(f, x)
        return logdensity(ν, y)
    end
end

function logdensity(pf::Pushforward{F}, y) where {F <: AbstractTransform}
    f = pf.f
    μ = pf.μ
    x = inverse(f, y)
    if pf.logjac
        _, ℓ = transform_and_logjac(f, x)
        return logdensity(μ, x) - ℓ
    else
        return logdensity(μ, x)
    end
end

Pullback(f::InverseTransform, ν, logjac=true) = Pushforward(f.transform, ν, logjac)

Pushforward(f::InverseTransform, ν, logjac=true) = Pullback(f.transform, ν, logjac)

Base.rand(rng::AbstractRNG, ν::Pushforward) = ν.f(rand(rng, ν.μ))

Base.rand(rng::AbstractRNG, μ::Pullback) = μ.f(rand(rng, μ.ν))

testvalue(ν::Pushforward) = transform(ν.f, testvalue(ν.μ))

testvalue(μ::Pullback) = transform(inverse(μ.f), testvalue(μ.ν))

basemeasure(μ::Pullback) = Pullback(μ.f, basemeasure(μ.ν), false)

basemeasure(ν::Pushforward) = Pushforward(ν.f, basemeasure(ν.μ), false)

TransformVariables.as(ν::Pushforward) = ν.f ∘ as(ν.μ)

TransformVariables.as(μ::Pullback) = inverse(μ.f) ∘ μ.ν

TransformVariables.as(::Lebesgue) = asℝ

# t = as𝕀
# μ = Normal()
# ν = Pushforward(t, μ)
# x = rand(μ)
# julia> logdensity(μ, x) ≈ logdensity(Pushforward(inverse(t), ν), x)
# true
