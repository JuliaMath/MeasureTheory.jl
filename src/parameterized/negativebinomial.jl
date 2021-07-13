# NegativeBinomial distribution

export NegativeBinomial
import Base
using StatsFuns
using SpecialFunctions

@parameterized NegativeBinomial(r,p) ≪ CountingMeasure(ℤ[0:∞])



(d::NegativeBinomial ≪ ::CountingMeasure{IntegerRange{a,b}}) where {a,b} = a ≤ 0 && b ≥ d.n

(::CountingMeasure{IntegerRange{a,b}} ≪ ::NegativeBinomial) where {a,b} = a ≥ 0 && b ≤ d.n

###############################################################################
@kwstruct NegativeBinomial(r, p)
    
function logdensity(d::NegativeBinomial{(:r, :p)}, y)
    (r, p) = (d.r, d.p)
    return -log(y + r) - logbeta(r, y+1) + y * log(p) + r * log1p(-p)
end

distproxy(d::NegativeBinomial{(:r, :p)}) = Dists.NegativeBinomial(d.r, d.p)

###############################################################################
@kwstruct NegativeBinomial(r, logitp)

function logdensity(d::NegativeBinomial{(:r, :logitp)}, y)
    (r, logitp) = (d.r, d.logitp)
    return -log(y + r) - logbeta(r, y+1) - y * log1pexp(-logitp) - r * log1pexp(logitp)
end

distproxy(d::NegativeBinomial{(:n,:logitp)}) = Dists.NegativeBinomial(d.n, logistic(d.logitp))

###############################################################################
@kwstruct NegativeBinomial(r, λ)
# mean λ, as in Poisson
# Converges to Poisson as r→∞

function logdensity(d::NegativeBinomial{(:r, :λ)}, y)
    (r, λ) = (d.r, d.λ)
    return -log(y + r) - logbeta(r, y+1) + y * log(λ) + r * log(r) - (y + r) * log(r + λ)
end

function distproxy(d::NegativeBinomial{(:r,:λ)})
    p = d.λ / (d.r + d.λ)
    return Dists.NegativeBinomial(d.r, p)
end

# # https://en.wikipedia.org/wiki/Negative_binomial_distribution#Gamma%E2%80%93Poisson_mixture
# function Base.rand(rng::AbstractRNG, d::NegativeBinomial{(:r,:λ)})
#     r = d.r
#     λ = d.λ
#     μ = rand(rng, Dists.Gamma(r, λ/r))
#     return rand(rng, Dists.Poisson(μ))
# end

###############################################################################
@kwstruct NegativeBinomial(r, logλ)

function logdensity(d::NegativeBinomial{(:r, :logλ)}, y)
    (r, logλ) = (d.r, d.logλ)
    λ = exp(logλ)
    return -log(y + r) - logbeta(r, y+1) + y * logλ + r * log(r) - (y + r) * log(r + λ)
end

function distproxy(d::NegativeBinomial{(:r,:logλ)})
    λ = exp(d.logλ)
    p = λ / (d.r + λ)
    return Dists.NegativeBinomial(d.r, p)
end

Base.rand(rng::AbstractRNG, d::NegativeBinomial{(:r,:logλ)}) = rand(rng, NegativeBinomial{(:r, :λ)}(d.r, exp(d.logλ)))

###############################################################################

asparams(::Type{<:NegativeBinomial}, ::Val{:p}) = as𝕀
asparams(::Type{<:NegativeBinomial}, ::Val{:logitp}) = asℝ
asparams(::Type{<:NegativeBinomial}, ::Val{:r}) = asℝ₊
asparams(::Type{<:NegativeBinomial}, ::Val{:λ}) = asℝ₊
