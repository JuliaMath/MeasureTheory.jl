# Beta distribution

export Beta

@parameterized Beta(α,β)

@kwstruct Beta(α, β)

@kwalias Beta [
    a     => α
    alpha => α
    b     => β
    beta  => β
]

TV.as(::Beta) = as𝕀

function logdensity(d::Beta{(:α, :β), Tuple{A,B}}, x::X) where {A,B,X}
    if static_hasmethod(xlogy, Tuple{A,X}) && static_hasmethod(xlog1py, Tuple{B,X})
        return xlogy(d.α - 1, x) + xlog1py(d.β - 1, -x) 
    else
        return (d.α - 1) * log(x) + (d.β - 1) * log1p(-x)
    end
end

function basemeasure(d::Beta{(:α,:β)})
    inbounds(x) = 0 < x < 1
    constℓ = 0.0
    varℓ() = - logbeta(d.α, d.β)
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Beta) = rand(rng, Dists.Beta(μ.α, μ.β))

distproxy(d::Beta{(:α, :β)}) = Dists.Beta(d.α, d.β)

asparams(::Type{<:Beta}, ::Val{:α}) = asℝ₊
asparams(::Type{<:Beta}, ::Val{:β}) = asℝ₊
