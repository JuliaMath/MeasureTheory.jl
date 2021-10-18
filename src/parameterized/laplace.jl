
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


function logdensity(d::Laplace{()} , x)
    return -abs(x)
end

Base.rand(rng::AbstractRNG, μ::Laplace{()}) = rand(rng, Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true

TV.as(::Laplace) = asℝ

distproxy(::Laplace{()}) = Dists.Laplace()
distproxy(d::Laplace{(:μ,)}) = Dists.Laplace(d.μ, 1.0)
distproxy(d::Laplace{(:σ,)}) = Dists.Laplace(0.0, d.σ)
distproxy(d::Laplace{(:μ,:σ)}) = Dists.Laplace(d.μ, d.σ)
distproxy(d::Laplace{(:ω,)}) = Dists.Laplace(0.0, inv(d.ω))
distproxy(d::Laplace{(:μ,:ω)}) = Dists.Laplace(d.μ, inv(d.ω))
