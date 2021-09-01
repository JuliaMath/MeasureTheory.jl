
# StudentT distribution

export StudentT

@parameterized StudentT(ν)

@kwstruct StudentT(ν)

StudentT(nt::NamedTuple{(:ν,:μ,:σ)}) = Affine(nt, StudentT())
StudentT(nt::NamedTuple{(:ν,:μ,:ω)}) = Affine(nt, StudentT())
StudentT(nt::NamedTuple{(:ν,:σ,)}) = Affine(nt, StudentT())
StudentT(nt::NamedTuple{(:ν,:ω,)}) = Affine(nt, StudentT())
StudentT(nt::NamedTuple{(:ν,:μ,)}) = Affine(nt, StudentT())

@affinepars StudentT

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

function logdensity(d::StudentT{(:ν,)}, x) 
    ν = d.ν
    return  - (ν + 1) / 2 * log1p(x^2 / ν)
end

function basemeasure(d::StudentT{(:ν,)})
    inbounds(x) = true
    constℓ = 0.0
    varℓ = loggamma((d.ν+1)/2) - loggamma(d.ν/2) - log(π * d.ν) / 2
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

TV.as(::StudentT) = asℝ

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

distproxy(d::StudentT{(:ν, :μ, :σ)}) = Dists.LocationScale(d.μ, d.σ, Dists.TDist(d.ν))

@half StudentT
@kwstruct StudentT()

HalfStudentT(ν,σ) = HalfStudentT(ν=ν,σ=σ)

asparams(::Type{<:StudentT}, ::Val{:ν}) = asℝ₊
