
# StudentT distribution

export StudentT, HalfStudentT

@parameterized StudentT(ν)

@kwstruct StudentT(ν)
@kwstruct StudentT(ν,μ)
@kwstruct StudentT(ν,σ)
@kwstruct StudentT(ν,μ,σ)

StudentT(nt::NamedTuple{(:ν,:μ,:σ)}) = affine(NamedTuple{(:μ,:σ)}(nt), StudentT(ν=nt.ν) ^ colsize(nt.σ))
StudentT(nt::NamedTuple{(:ν,:μ,:ω)}) = affine(NamedTuple{(:μ,:ω)}(nt), StudentT(ν=nt.ν) ^ rowsize(nt.ω))
StudentT(nt::NamedTuple{(:ν,:σ,)}) = affine(NamedTuple{(:σ,)}(nt), StudentT(ν=nt.ν) ^ colsize(nt.σ))
StudentT(nt::NamedTuple{(:ν,:ω,)}) = affine(NamedTuple{(:ω,)}(nt), StudentT(ν=nt.ν) ^ rowsize(nt.ω))
StudentT(nt::NamedTuple{(:ν,:μ,)}) = affine(NamedTuple{(:μ,)}(nt), StudentT(ν=nt.ν) ^ size(nt.μ))

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

function logdensity(d::StudentT{(:ν,), Tuple{T}}, x::Number) where {T<:Number}
    ν = d.ν
    return  xlog1py((ν + 1) / (-2), x^2 / ν)
end

function logdensity(d::StudentT{(:ν,)}, x) 
    ν = d.ν
    return  (ν + 1) / (-2) * log1p(x^2 / ν)
end

function basemeasure(d::StudentT{(:ν,)})
    inbounds(x) = true
    constℓ = 0.0
    varℓ() = loggamma((d.ν+1)/2) - loggamma(d.ν/2) - log(π * d.ν) / 2
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

TV.as(::StudentT) = asℝ

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

distproxy(d::StudentT{(:ν,)}) = Dists.TDist(d.ν)

@half StudentT

HalfStudentT(ν, σ) = HalfStudentT((ν=ν, σ=σ))

asparams(::Type{<:StudentT}, ::Val{:ν}) = asℝ₊
