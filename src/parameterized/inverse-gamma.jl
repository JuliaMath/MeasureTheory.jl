
# InverseGamma distribution

using SpecialFunctions
using StatsFuns
export InverseGamma

@parameterized InverseGamma(α) ≃ Lebesgue(ℝ₊)

function logdensity(μ::InverseGamma{(:α,)}, x) 
    α = μ.α
    xinv = 1/x

    return (α + 1) * log(xinv) - xinv - loggamma(α) 
end

Base.rand(rng::AbstractRNG, T::Type, μ::InverseGamma{(:shape,)}) = rand(rng, Dists.InverseGamma(μ.shape))

≪(::InverseGamma, ::Lebesgue{X}) where X <: Real = true

TV.as(::InverseGamma) = asℝ₊

asparams(::Type{<:InverseGamma}, ::Val{:α}) = asℝ₊

@μσ_methods InverseGamma(α)
