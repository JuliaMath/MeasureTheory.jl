
# StudentT distribution

export StudentT, HalfStudentT

@parameterized StudentT(ν)

@kwstruct StudentT(ν)
@kwstruct StudentT(ν, μ)
@kwstruct StudentT(ν, σ)
@kwstruct StudentT(ν, ω)
@kwstruct StudentT(ν, μ, σ)
@kwstruct StudentT(ν, μ, ω)

for N in AFFINEPARS
    @eval begin
        proxy(d::StudentT{(:ν, $N...)}) =
            affine(NamedTuple{$N}(params(d)), StudentT((ν = d.ν,)))
    end
end

logdensity_def(d::StudentT, x) = logdensity_def(proxy(d), x)
basemeasure(d::StudentT) = basemeasure(proxy(d))
Base.rand(rng::AbstractRNG, ::Type{T}, μ::StudentT) where {T} = rand(rng, T, proxy(μ))

# @affinepars StudentT

StudentT(ν, μ, σ) = StudentT((ν = ν, μ = μ, σ = σ))

@kwalias StudentT [
    df       => ν
    nu       => ν
    location => μ
    mean     => μ # Doesn't really always exist, but ok fine
    mu       => μ
    scale    => σ
    sigma    => σ
]

@inline function logdensity_def(d::StudentT{(:ν,),Tuple{T}}, x::Number) where {T<:Number}
    ν = d.ν
    return xlog1py((ν + 1) / (-2), x^2 / ν)
end

@inline function logdensity_def(d::StudentT{(:ν,)}, x)
    ν = d.ν
    return (ν + 1) / (-2) * log1p(x^2 / ν)
end

@inline function basemeasure(d::StudentT{(:ν,)})
    constℓ = 0.0
    varℓ() = loggamma((d.ν + 1) / 2) - loggamma(d.ν / 2) - log(π * d.ν) / 2
    base = Lebesgue(ℝ)
    FactoredBase(constℓ, varℓ, base)
end

xform(::StudentT) = asℝ

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

proxy(d::StudentT{(:ν,)}) = Dists.TDist(d.ν)

@half StudentT

HalfStudentT(ν, σ) = HalfStudentT((ν = ν, σ = σ))

asparams(::Type{<:StudentT}, ::StaticSymbol{:ν}) = asℝ₊

insupport(::StudentT, x) = true
insupport(::StudentT) = Returns(true)
