"""
    struct ScaledMeasure{R,M} <: AbstractMeasure
        logscale :: R
        base :: M
    end
    

"""
struct ScaledMeasure{R,M} <: AbstractMeasure
    logscale :: R
    base :: M
end

function logdensity(sm::ScaledMeasure{R,M}, x::X) where {X, R, M <: AbstractMeasure}
    logdensity(sm.base, x) + sm.logscale
end

function Base.:*(k::Number, m::AbstractMeasure) 
    if iszero(k)
        return TrivialMeasure
    else
        return ScaledMeasure{typeof(k), typeof(m)}(log(k),m)
    end
end

Base.:*(m::AbstractMeasure, k::Real) = k * m

≪(::M, ::ScaledMeasure{R,M}) where {R,M,X} = true
≪(::ScaledMeasure{R,M}, ::M) where {R,M,X} = true

baseMeasure(μ::ScaledMeasure{R,M}) where {R,M,X} = μ.base

sampletype(μ:: ScaledMeasure) = sampletype(μ.base)
