
# StudentT distribution

export StudentT, HalfStudentT

@parameterized StudentT(ν)

@kwstruct StudentT(ν)
@kwstruct StudentT(ν,μ)
@kwstruct StudentT(ν,σ)
@kwstruct StudentT(ν,ω)
@kwstruct StudentT(ν,μ,σ)
@kwstruct StudentT(ν,μ,ω)

for N in AFFINEPARS
    @eval begin
        proxy(d::StudentT{(:ν,$N...)}) = affine(NamedTupleTools.select(params(d), $N), StudentT((ν=d.ν,)))
        logdensity(d::StudentT{(:ν, $N...)}, x) = logdensity(proxy(d), x)
        basemeasure(d::StudentT{(:ν, $N...)}) = basemeasure(proxy(d))
    end
end

# @affinepars StudentT

StudentT(ν, μ, σ) = StudentT((ν=ν, μ=μ, σ=σ))

@kwalias StudentT [
    df       => ν
    nu       => ν
    location => μ
    mean     => μ # Doesn't really always exist, but ok fine
    mu       => μ
    scale    => σ
    sigma    => σ
]

@inline function logdensity(d::StudentT{(:ν,), Tuple{T}}, x::Number) where {T<:Number}
    ν = d.ν
    return  xlog1py((ν + 1) / (-2), x^2 / ν)
end

@inline function logdensity(d::StudentT{(:ν,)}, x) 
    ν = d.ν
    return  (ν + 1) / (-2) * log1p(x^2 / ν)
end

@inline function basemeasure(d::StudentT{(:ν,)})
    inbounds = Returns(true)
    constℓ = 0.0
    varℓ() = loggamma((d.ν+1)/2) - loggamma(d.ν/2) - log(π * d.ν) / 2
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

TV.as(::StudentT) = asℝ

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

distproxy(d::StudentT{(:ν,)}) = Dists.TDist(d.ν)
distproxy(d::StudentT{(:ν,:μ)}) = Dists.LocationScale(d.μ, 1.0, Dists.TDist(d.ν))
distproxy(d::StudentT{(:ν,:σ)}) = Dists.LocationScale(0.0, d.σ, Dists.TDist(d.ν))
distproxy(d::StudentT{(:ν,:ω)}) = Dists.LocationScale(0.0, inv(d.ω), Dists.TDist(d.ν))
distproxy(d::StudentT{(:ν,:μ,:σ)}) = Dists.LocationScale(d.μ, d.σ, Dists.TDist(d.ν))
distproxy(d::StudentT{(:ν,:μ,:ω)}) = Dists.LocationScale(d.μ, inv(d.ω), Dists.TDist(d.ν))


@half StudentT

HalfStudentT(ν, σ) = HalfStudentT((ν=ν, σ=σ))

asparams(::Type{<:StudentT}, ::Val{:ν}) = asℝ₊
