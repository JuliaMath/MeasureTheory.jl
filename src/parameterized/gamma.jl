# Gamma distribution

export Gamma

@parameterized Gamma(α,β) ≪ Lebesgue(ℝ₊)

@kwstruct Gamma(α, β)

@kwalias Gamma [
    a     => α
    alpha => α
    b     => β
    beta  => β
]

TV.as(::Gamma) = asℝ₊

function logdensity(d::Gamma{(:α, :β)}, x)
    return d.α * log(d.β) + (d.α - 1) * log(x) - d.β * x - loggamma(d.α)
end

distproxy(d::Gamma{(:α, :β)}) = Dists.Gamma(d.α, 1/d.β)

asparams(::Type{<:Gamma}, ::Val{:α}) = asℝ₊
asparams(::Type{<:Gamma}, ::Val{:β}) = asℝ₊
