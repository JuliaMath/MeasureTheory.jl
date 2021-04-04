# NegativeBinomial distribution

export NegativeBinomial
import Base
using StatsFuns
using SpecialFunctions

@measure NegativeBinomial(r,p)

basemeasure(μ::NegativeBinomial) = CountingMeasure(ℤ[0:∞])

(d::NegativeBinomial ≪ ::CountingMeasure{IntegerRange{a,b}}) where {a,b} = a ≤ 0 && b ≥ d.n

(::CountingMeasure{IntegerRange{a,b}} ≪ ::NegativeBinomial) where {a,b} = a ≥ 0 && b ≤ d.n

###############################################################################
# (r, p)
    
function logdensity(d::NegativeBinomial{(:r, :p)}, y)
    (r, p) = (d.r, d.p)
    return -log(y + r) - logbeta(r, y+1) + y * log(p) + r * log1p(-p)
end

###############################################################################
# (r, logitp)

function logdensity(d::NegativeBinomial{(:r, :logitp)}, y)
    (r, logitp) = (d.r, d.logitp)
    return -log(y + r) - logbeta(r, y+1) - y * log1pexp(-logitp) - r * log1pexp(logitp)
end

###############################################################################
# (r,λ) 
# mean λ, as in Poisson
# Converges to Poisson as r→∞

function logdensity(d::NegativeBinomial{(:r, :λ)}, y)
    (r, λ) = (d.r, d.λ)
    return -log(y + r) - logbeta(r, y+1) + y * log(λ) + r * log(r) - (y + r) * log(r + λ)
end

representative(d::NegativeBinomial) = CountingMeasure(ℤ[0:d.n])

distproxy(d::NegativeBinomial{(:r, :p)}) = Dists.NegativeBinomial(d.r, d.p)
distproxy(d::NegativeBinomial{(:n,:logitp)}) = Dists.NegativeBinomial(d.n, logistic(d.logitp))

function distproxy(d::NegativeBinomial{(:r,:λ)})
    p = d.λ / (d.r + d.λ)
    return Dists.NegativeBinomial(d.r, p)
end

# https://en.wikipedia.org/wiki/Negative_binomial_distribution#Gamma%E2%80%93Poisson_mixture
function Base.rand(rng::AbstractRNG, d::NegativeBinomial{(:r,:λ)})
    r = d.r
    λ = d.λ
    μ = rand(rng, Dists.Gamma(r, λ/r))
    return rand(rng, Dists.Poisson(μ))
end

asparams(::Type{<:NegativeBinomial}, ::Val{:p}) = as𝕀
asparams(::Type{<:NegativeBinomial}, ::Val{:logitp}) = asℝ
asparams(::Type{<:NegativeBinomial}, ::Val{:r}) = asℝ₊
asparams(::Type{<:NegativeBinomial}, ::Val{:λ}) = asℝ₊
