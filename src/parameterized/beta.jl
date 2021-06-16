# Beta distribution

import StatsFuns
export Beta

@parameterized Beta(α,β) ≃ Lebesgue(𝕀)

@kwstruct Beta(α, β)

@kwalias Beta [
    a     => α
    alpha => α
    b     => β
    beta  => β
]

TV.as(::Beta) = as𝕀

function logdensity(d::Beta{(:α, :β)}, x)
    return (d.α - 1) * log(x) + (d.β - 1) * log(1 - x) - logbeta(d.α, d.β)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Beta) = rand(rng, Dists.Beta(μ.α, μ.β))

distproxy(d::Beta{(:α, :β)}) = Dists.Beta(d.α, d.β)
