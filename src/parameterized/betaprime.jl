# BetaPrime distribution
# Some code based on https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/continuous/betaprime.jl

export BetaPrime

import StatsFuns

@parameterized BetaPrime(α, β)

@kwstruct BetaPrime(α, β)

@kwalias BetaPrime [
    a => α
    alpha => α
    b => β
    betaPrime => β
]

function logdensity_def(d::BetaPrime{(:α, :β)}, x::Real)
    α = d.α
    β = d.β
    _x = max(0, x)
    z = xlogy(α - 1, _x) - (α + β) * log1p(_x) - logbeta(α, β)
    return x < 0 ? oftype(z, -Inf) : z
end


function zval(::Type{<:BetaPrime}, x::Real)
    y = max(x, 0)
    z = y / (1 + y)
    # map `Inf` to `Inf` (otherwise it returns `NaN`)
    return isinf(x) && x > 0 ? oftype(z, Inf) : z
end

function xval(::Type{<:BetaPrime}, z::Real)
    return z / (1 - z)
end


MeasureBase.basemeasure(d::BetaPrime{(:α, :β)}) = LebesgueBase()

function Base.rand(rng::AbstractRNG, ::Type{T}, μ::BetaPrime) where {T}
    rand(rng, T, Gamma(μ.α)) / rand(rng, T, Gamma(μ.β))
end

insupport(::BetaPrime, x) = x > 0

function smf(d::BetaPrime{(:α, :β)}, x::Real)
    smf(Beta(d.α, d.β), zval(BetaPrime, x))
end

function invsmf(d::BetaPrime{(:α, :β)}, p)
    xval(BetaPrime, invsmf(Beta(d.α, d.β), p))
end

proxy(d::BetaPrime{(:α, :β)}) = Dists.BetaPrime(d.α, d.β)
