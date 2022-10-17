# Dirichlet distribution
using MeasureBase: Simplex

export Dirichlet

using FillArrays

@parameterized Dirichlet(α)

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
    mapreduce(+, d.α, x) do αj, xj
        xlogy(αj - 1, xj)
    end
end

Base.rand(rng::AbstractRNG, T::Type, μ::Dirichlet) = rand(rng, Dists.Dirichlet(μ.α))


function testvalue(::Type{T}, d::Dirichlet{(:α,)}) where {T}
    n = length(d.α)
    Fill(inv(convert(T, n)), n)
end

@inline function insupport(d::Dirichlet{(:α,)}, x)
    length(x) == length(d.α) && sum(x) ≈ 1.0
end

as(μ::Dirichlet) = TV.UnitSimplex(length(μ.α))
