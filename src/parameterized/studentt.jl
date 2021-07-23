
# StudentT distribution

using SpecialFunctions
using StatsFuns
export StudentT

@parameterized StudentT(ν) ≪ (1/sqrtπ) * Lebesgue(ℝ)

@kwstruct StudentT(ν)

@μσ_methods StudentT(ν)

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
    halfνp1 = (ν+1)/2
    return loggamma(halfνp1) - loggamma(ν/2) + ν * log(ν) - halfνp1 * log(x^2 + ν)
end

TV.as(::StudentT) = asℝ

Base.rand(rng::AbstractRNG, T::Type, μ::StudentT{(:ν,)}) = rand(rng, Dists.TDist(μ.ν))

distproxy(d::StudentT{(:ν, :μ, :σ)}) = Dists.LocationScale(d.μ, d.σ, Dists.TDist(d.ν))

@half StudentT(ν)
@kwstruct StudentT()

asparams(::Type{<:StudentT}, ::Val{:ν}) = asℝ₊
