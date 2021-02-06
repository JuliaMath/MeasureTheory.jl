# Beta distribution

import StatsFuns
export Beta

@measure Beta(α,β) ≃ Lebesgue(𝕀)

function logdensity(d::Beta{(:α, :β)}, x)
    return (d.α - 1) * log(x) + (d.β - 1) * log(1 - x) - logbeta(d.α, d.β)
end

Base.rand(rng::AbstractRNG, μ::Beta) = rand(rng, Dists.Beta(μ.α, μ.β))

≪(::Beta, ::Lebesgue{X}) where X <: Real = true
representative(::Beta) = Lebesgue(𝕀)
