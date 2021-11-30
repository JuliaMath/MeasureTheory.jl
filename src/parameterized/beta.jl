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

TV.as(::Beta) = as𝕀

@inline function logdensity_def(d::Beta{(:α, :β),Tuple{A,B}}, x::X) where {A,B,X}
    return xlogy(d.α - 1, x) + xlog1py(d.β - 1, -x)
end

@inline function basemeasure(d::Beta{(:α, :β)})
    inbounds(x) = 0 < x < 1
    constℓ = 0.0
    varℓ() = -logbeta(d.α, d.β)
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

basemeasure_depth(::Beta) = static(2)
basemeasure_depth(::Type{T}) where {T<:Beta} = static(2)

Base.rand(rng::AbstractRNG, T::Type, μ::Beta) = rand(rng, Dists.Beta(μ.α, μ.β))

distproxy(d::Beta{(:α, :β)}) = Dists.Beta(d.α, d.β)

asparams(::Type{<:Beta}, ::StaticSymbol{:α}) = asℝ₊
asparams(::Type{<:Beta}, ::StaticSymbol{:β}) = asℝ₊
