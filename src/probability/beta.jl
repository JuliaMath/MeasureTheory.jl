
# Beta distribution

import StatsFuns
export Beta


@measure Beta(α,β) ≃ Lebesgue(𝕀)


# Standard Beta

function logdensity(d::Beta{(:α, :β)}, x::MaybeSym{T}) where {T <: Number}
    return (d.α - 1) * log(x) + (d.β - 1) * log(1 - x) - logbeta(d.α, d.β)
end

# Beta() = Beta{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::Beta) = Real

Base.rand(rng::AbstractRNG, μ::Beta) = rand(rng, Dists.Beta(μ.α, μ.β))

≪(::Beta, ::Lebesgue{X}) where X <: Real = true
representative(::Beta) = Lebesgue(𝕀)
