
# Laplace distribution

export Laplace

@parameterized Laplace() ≪ (1/2) * Lebesgue(ℝ)

Laplace(nt::NamedTuple{(:μ,:σ)}) = Affine(nt, Laplace())
Laplace(nt::NamedTuple{(:μ,:ω)}) = Affine(nt, Laplace())
Laplace(nt::NamedTuple{(:σ,)}) = Affine(nt, Laplace())
Laplace(nt::NamedTuple{(:ω,)}) = Affine(nt, Laplace())
Laplace(nt::NamedTuple{(:μ,)}) = Affine(nt, Laplace())

@affinepars Laplace


function logdensity(d::Laplace{()} , x)
    return -abs(x)
end

Base.rand(rng::AbstractRNG, μ::Laplace{()}) = rand(rng, Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true

TV.as(::Laplace) = asℝ

# @μσ_methods Laplace()
@half Laplace

HalfLaplace(σ) = HalfLaplace(σ=σ)

distproxy(::Laplace{()}) = Dists.Laplace()
