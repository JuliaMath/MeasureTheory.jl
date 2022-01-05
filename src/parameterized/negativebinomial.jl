# NegativeBinomial distribution

export NegativeBinomial
import Base

@parameterized NegativeBinomial(r, p)

basemeasure(::NegativeBinomial) = Counting(BoundedInts(static(0), static(Inf)))

function Base.rand(rng::AbstractRNG, ::Type{T}, d::NegativeBinomial) where {T}
    rand(rng, proxy(d))
end

###############################################################################
@kwstruct NegativeBinomial(r, p)

@inline function logdensity_def(d::NegativeBinomial{(:r, :p)}, y)
    (r, p) = (d.r, d.p)
    return -log(y + r) - logbeta(r, y + 1) + xlogy(y, p) + xlog1py(r, -p)
end

proxy(d::NegativeBinomial{(:r, :p)}) = Dists.NegativeBinomial(d.r, d.p)

###############################################################################
@kwstruct NegativeBinomial(r, logitp)

@inline function logdensity_def(d::NegativeBinomial{(:r, :logitp)}, y)
    (r, logitp) = (d.r, d.logitp)
    return -log(y + r) - logbeta(r, y + 1) - y * log1pexp(-logitp) - r * log1pexp(logitp)
end

proxy(d::NegativeBinomial{(:n, :logitp)}) = Dists.NegativeBinomial(d.n, logistic(d.logitp))

###############################################################################
@kwstruct NegativeBinomial(r, λ)
# mean λ, as in Poisson
# Converges to Poisson as r→∞

@inline function logdensity_def(d::NegativeBinomial{(:r, :λ)}, y)
    (r, λ) = (d.r, d.λ)
    return -log(y + r) - logbeta(r, y + 1) + xlogy(y, λ) + xlogx(r) - xlogy(y + r, r + λ)
end

function proxy(d::NegativeBinomial{(:r, :λ)})
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

@inline function logdensity_def(d::NegativeBinomial{(:r, :logλ)}, y)
    (r, logλ) = (d.r, d.logλ)
    λ = exp(logλ)
    return -log(y + r) - logbeta(r, y + 1) + y * logλ + xlogx(r) - xlogy(y + r, r + λ)
end

function proxy(d::NegativeBinomial{(:r, :logλ)})
    λ = exp(d.logλ)
    p = λ / (d.r + λ)
    return Dists.NegativeBinomial(d.r, p)
end

###############################################################################

asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:p}) = as𝕀
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:logitp}) = asℝ
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:r}) = asℝ₊
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:λ}) = asℝ₊
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:logλ}) = asℝ
