
# Laplace distribution

export Laplace

@parameterized Laplace() ≪ (1/2) * Lebesgue(ℝ)

Laplace(nt::NamedTuple{(:μ,:σ)}) = affine(nt, Laplace())
Laplace(nt::NamedTuple{(:μ,:ω)}) = affine(nt, Laplace())
Laplace(nt::NamedTuple{(:σ,)}) = affine(nt, Laplace())
Laplace(nt::NamedTuple{(:ω,)}) = affine(nt, Laplace())
Laplace(nt::NamedTuple{(:μ,)}) = affine(nt, Laplace())

@affinepars Laplace


function logdensity(d::Laplace{()} , x)
    return -abs(x)
end

Base.rand(rng::AbstractRNG, μ::Laplace{()}) = rand(rng, Dists.Laplace())

≪(::Laplace, ::Lebesgue{X}) where X <: Real = true

TV.as(::Laplace) = asℝ

distproxy(::Laplace{()}) = Dists.Laplace()
