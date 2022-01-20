
# Laplace distribution

export Laplace

@parameterized Laplace() 

@kwstruct Laplace()
@kwstruct Laplace(μ)
@kwstruct Laplace(σ)
@kwstruct Laplace(μ, σ)
@kwstruct Laplace(ω)
@kwstruct Laplace(μ, ω)

for N in AFFINEPARS
    @eval begin
        proxy(d::Laplace{$N}) = affine(params(d), Laplace())
        @useproxy Laplace{$N}
    end
end



@inline function logdensity_def(d::Laplace{()}, x)
    return -abs(x)
end

logdensity_def(d::Laplace, x) = logdensity_def(proxy(d), x)


basemeasure(::Laplace{()}) = WeightedMeasure(static(-logtwo), Lebesgue(ℝ))

# @affinepars Laplace


Base.rand(rng::AbstractRNG, ::Type{T}, μ::Laplace{()}) where {T} =
    rand(rng, Dists.Laplace())
Base.rand(rng::AbstractRNG, ::Type{T}, μ::Laplace) where {T} = Base.rand(rng, T, proxy(μ))

≪(::Laplace, ::Lebesgue{X}) where {X<:Real} = true

TV.as(::Laplace) = asℝ
