export ScaledInvChiSq
using MeasureTheory: @parameterized
using KeywordCalls: @kwstruct
using MeasureBase: params

@parameterized ScaledInvChiSq(ν, s²)

@kwstruct ScaledInvChiSq(ν, s²)

MeasureBase.basemeasure(::ScaledInvChiSq) = LebesgueBase()

function MeasureBase.insupport(d::ScaledInvChiSq, x)
    x > 0
end

function MeasureBase.proxy(d::ScaledInvChiSq{(:ν, :s²)})
    pars = params(d)
    halfν = 0.5 * pars.ν
    s² = pars.s²
    α = halfν
    β = halfν * s²
    return Dists.InverseGamma(α, β)
end

MeasureBase.logdensity_def(d::ScaledInvChiSq{(:ν, :s²)}, x) = logdensityof(proxy(d), x)

function Base.rand(rng::AbstractRNG, ::Type{T}, d::ScaledInvChiSq) where {T}
    return rand(rng, proxy(d))
end
