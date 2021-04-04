
# Normal distribution

import StatsFuns
export Normal

@measure Normal(μ,σ)

# To pass `@code_warntype`
const _Normal_basemeasure = (1/sqrt2π) * Lebesgue(ℝ)

basemeasure(::Normal)=  _Normal_basemeasure

logdensity(d::Normal{()} , x) = - x^2 / 2 

Base.rand(rng::Random.AbstractRNG, T::Type, μ::Normal{()}) = randn(rng, T)

@μσ_methods Normal()

@half Normal()

HalfNormal(σ) = HalfNormal(σ = σ)


###############################################################################
# (μ,σ²)

function logdensity(d::Normal{(:σ²)}, x)
    σ² = d.σ²
    -0.5 * (log(σ²) + (x^2/σ²))
end

function logdensity(d::Normal{(:μ,:σ²)}, x)
    μ = d.μ
    σ² = d.σ²
    -0.5 * (log(σ²) + ((x - μ)^2/σ²))
end

###############################################################################
# (μ,τ)

function logdensity(d::Normal{(:τ)}, x)
    τ = d.τ
    0.5 * (log(τ) - τ * x^2)
end

function logdensity(d::Normal{(:μ,:τ)}, x)
    μ = d.μ
    τ = d.τ
    0.5 * (log(τ) - τ * (x - μ)^2)
end


###############################################################################
# (μ,logσ)

function logdensity(d::Normal{(:μ,:logσ)}, x)
    μ = d.μ
    logσ = d.logσ
    -logσ - 0.5(exp(-2logσ)*((x - μ)^2))
end



distproxy(d::Normal{(:μ, :σ)}) = Dists.Normal(d.μ, d.σ)

asparams(::Type{<:Normal}, ::Val{:μ}) = asℝ
asparams(::Type{<:Normal}, ::Val{:σ}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:logσ}) = asℝ
asparams(::Type{<:Normal}, ::Val{:σ²}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:τ}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:logτ}) = asℝ
