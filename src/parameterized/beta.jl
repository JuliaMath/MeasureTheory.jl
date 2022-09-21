# Beta distribution

export Beta

@parameterized Beta(α, β)

@kwstruct Beta(α, β)

@kwalias Beta [
    a => α
    alpha => α
    b => β
    beta => β
]

@inline function logdensity_def(d::Beta{(:α, :β),Tuple{A,B}}, x::X) where {A,B,X}
    return xlogy(d.α - 1, x) + xlog1py(d.β - 1, -x)
end

@inline function basemeasure(d::Beta{(:α, :β)})
    ℓ = -logbeta(d.α, d.β)
    weightedmeasure(ℓ, LebesgueBase())
end

Base.rand(rng::AbstractRNG, T::Type, μ::Beta) = rand(rng, Dists.Beta(μ.α, μ.β))

proxy(d::Beta{(:α, :β)}) = Dists.Beta(d.α, d.β)

insupport(::Beta, x) = in𝕀(x)
insupport(::Beta) = in𝕀
