
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
    return - (ν+1)/2 * log(x^2 + ν)
end

function basemeasure(d::StudentT{(:ν,)})
    function ℓ(nt)
        ν = nt.ν
        loggamma((ν+1)/2) - loggamma(ν/2) + ν * log(ν)
    end
    return ParamWeightedMeasure(ℓ, params(d), (1/sqrtπ) * Lebesgue(ℝ))
end

# TV.as(::StudentT) = asℝ

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

distproxy(d::StudentT{(:ν, :μ, :σ)}) = Dists.LocationScale(d.μ, d.σ, Dists.TDist(d.ν))

@half StudentT
@kwstruct StudentT()

HalfStudentT(ν,σ) = HalfStudentT(ν=ν,σ=σ)

asparams(::Type{<:StudentT}, ::Val{:ν}) = asℝ₊
