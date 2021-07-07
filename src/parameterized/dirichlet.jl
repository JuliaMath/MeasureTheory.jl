# Dirichlet distribution

import StatsFuns
export Dirichlet

@parameterized Dirichlet(α)

TV.as(d::Dirichlet{(:α,)}) = TV.UnitSimplex(length(d.α))

function basemeasure(μ::Dirichlet{(:α,)})
    t = as(μ)
    d = dimension(t)
    return Pushforward(t, Lebesgue(ℝ)^d, false)
end

@kwstruct Dirichlet(α)

function logdensity(d::Dirichlet{(:α,)}, x)
    α = d.α
    s = 0.0
    for j in eachindex(x)
        s += α[j] * log(x[j])
    end
    return s
end

Base.rand(rng::AbstractRNG, T::Type, μ::Dirichlet) = rand(rng, Dists.Dirichlet(μ.α))

distproxy(d::Dirichlet{(:α,)}) = Dists.Dirichlet(d.α)
