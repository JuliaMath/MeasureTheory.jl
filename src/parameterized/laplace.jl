
# Laplace distribution

export Laplace

@parameterized Laplace(μ,σ) ≪ (1/2) * Lebesgue(ℝ)

@kwstruct Laplace()

function logdensity(d::Laplace{()} , x)
    return -abs(x)
end

Base.rand(rng::AbstractRNG, μ::Laplace{()}) = rand(rng, Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true

# TV.as(::Laplace) = asℝ

@μσ_methods Laplace()
@half Laplace()
@kwstruct HalfLaplace()

@σ_methods HalfLaplace()
HalfLaplace(σ) = HalfLaplace(σ=σ)

distproxy(::Laplace{()}) = Dists.Laplace()
