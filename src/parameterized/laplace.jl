
# Laplace distribution

export Laplace

@parameterized Laplace() ≪ (1/2) * Lebesgue(ℝ)

for N in AFFINEPARS
    @eval begin
        proxy(d::Laplace{$N}) = affine(params(d), Laplace())
        logdensity(d::Laplace{$N}, x) = logdensity(proxy(d), x)
        basemeasure(d::Laplace{$N}) = basemeasure(proxy(d))
    end
end

# @affinepars Laplace


@inline function logdensity(d::Laplace{()} , x)
    return -abs(x)
end

Base.rand(rng::AbstractRNG, ::Type{T}, μ::Laplace{()}) where {T} = rand(rng, Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true

TV.as(::Laplace) = asℝ
