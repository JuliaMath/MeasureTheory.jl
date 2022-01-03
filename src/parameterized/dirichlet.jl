# Dirichlet distribution
using MeasureBase: Simplex

export Dirichlet

using FillArrays

@parameterized Dirichlet(α)

TV.as(d::Dirichlet{(:α,)}) = TV.UnitSimplex(length(d.α))

@inline function basemeasure(μ::Dirichlet{(:α,)})
    α = μ.α
    t = as(μ)
    d = TV.dimension(t)
    logw = loggamma(sum(α)) - sum(loggamma, α)
    return WeightedMeasure(logw, Lebesgue(Simplex()))
end

@kwstruct Dirichlet(α)

Dirichlet(k::Integer, α) = Dirichlet(Fill(α, k))

@inline function logdensity_def(d::Dirichlet{(:α,)}, x)
    α = d.α
    s = 0.0
    for j in eachindex(x)
        s += xlogy(α[j] - 1, x[j])
    end
    return s
end

Base.rand(rng::AbstractRNG, T::Type, μ::Dirichlet) = rand(rng, Dists.Dirichlet(μ.α))

distproxy(d::Dirichlet{(:α,)}) = Dists.Dirichlet(d.α)

function testvalue(d::Dirichlet{(:α,)})
    n = length(d.α)
    Fill(1 / n, n)
end
