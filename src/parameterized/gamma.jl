# Gamma distribution

# REFERENCES
# https://mc-stan.org/docs/2_27/functions-reference/gamma-distribution.html
# https://en.wikipedia.org/wiki/Gamma_distribution
# https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Gamma

export Gamma

@parameterized Gamma(α, θ) ≪ Lebesgue(ℝ₊)

@kwstruct Gamma(α, θ)
@kwstruct Gamma(α, β)

@kwalias Gamma [
    a     => α
    alpha => α
    shape => α
    b     => β
    beta  => β
    rate  => β
    theta => θ
    scale => θ
]

TV.as(::Gamma) = asℝ₊

function logdensity(d::Gamma{(:α, :θ)}, x)
    return (d.α - 1) * log(x) - x / d.θ - d.α * log(d.θ) - loggamma(d.α)
end

function logdensity(d::Gamma{(:α, :β)}, x)
    return (d.α - 1) * log(x) - x * d.β + d.α * log(d.β) - loggamma(d.α)
end

distproxy(d::Gamma{(:α, :θ)}) = Dists.Gamma(d.α, d.θ)
distproxy(d::Gamma{(:α, :β)}) = Dists.Gamma(d.α, 1/d.β)

asparams(::Type{<:Gamma}, ::Val{:α}) = asℝ₊
asparams(::Type{<:Gamma}, ::Val{:β}) = asℝ₊
asparams(::Type{<:Gamma}, ::Val{:θ}) = asℝ₊
