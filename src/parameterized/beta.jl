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

function tbasemeasure_type(::Type{Beta{(:α, :β), T}}) where {T}
    R = float(eltype(T))
    FactoredBase{Base.Fix2{typeof(in), MeasureBase.BoundedReals{Static.StaticFloat64{0.0}, Static.StaticFloat64{1.0}}}, R, Returns{R}, Lebesgue{MeasureBase.RealNumbers}}
end

@inline function basemeasure(d::Beta{(:α, :β)})
    inbounds = in(𝕀)
    constℓ = 0.0
    varℓ = Returns(-logbeta(d.α, d.β))
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

Base.rand(rng::AbstractRNG, T::Type, μ::Beta) = rand(rng, Dists.Beta(μ.α, μ.β))

distproxy(d::Beta{(:α, :β)}) = Dists.Beta(d.α, d.β)

asparams(::Type{<:Beta}, ::StaticSymbol{:α}) = asℝ₊
asparams(::Type{<:Beta}, ::StaticSymbol{:β}) = asℝ₊
