
# Laplace distribution

export Laplace

@parameterized Laplace()

@kwstruct Laplace()
@kwstruct Laplace(μ)
@kwstruct Laplace(σ)
@kwstruct Laplace(μ, σ)
@kwstruct Laplace(λ)
@kwstruct Laplace(μ, λ)

for N in AFFINEPARS
    @eval begin
        proxy(d::Laplace{$N}) = affine(params(d), Laplace())
        @useproxy Laplace{$N}
    end
end

insupport(::Laplace, x) = true

@inline function logdensity_def(d::Laplace{()}, x)
    return -abs(x)
end

Laplace(μ, σ) = Laplace((μ=μ, σ=σ))

basemeasure(::Laplace{()}) = WeightedMeasure(static(-logtwo), LebesgueBase())

# @affinepars Laplace

function Base.rand(rng::AbstractRNG, ::Type{T}, μ::Laplace{()}) where {T}
    sign = rand(rng, Bool)
    absx = randexp(rng, T)
    sign == true ? absx : -absx
end

Base.rand(rng::AbstractRNG, ::Type{T}, μ::Laplace) where {T} = Base.rand(rng, T, proxy(μ))

≪(::Laplace, ::Lebesgue{X}) where {X<:Real} = true
