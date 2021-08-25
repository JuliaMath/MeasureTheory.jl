# Gumbel distribution

export Gumbel

@parameterized Gumbel() ≃ Lebesgue(ℝ)

@kwstruct Gumbel()

Gumbel(nt::NamedTuple{(:σ,)}) = Affine(nt, Gumbel())

@affinepars Gumbel


function logdensity(d::Gumbel{()} , x)
    return -exp(-x) - x
end

import Base

function Base.rand(rng::AbstractRNG, d::Gumbel{()})
    u = rand(rng)
    -log(-log(u))
end

# TV.as(::Gumbel) = asℝ

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true

distproxy(::Gumbel{()}) = Dists.Gumbel()
