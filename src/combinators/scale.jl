"""
    struct ScaledMeasure{R,M,X} <: AbstractMeasure{X}
        logscale :: R
        base :: M
    end
    

"""
struct ScaledMeasure{R,M,X} <: AbstractMeasure{X}
    logscale :: R
    base :: M
end

function logdensity(sm::ScaledMeasure{R,M,X}, x::X) where {X, R <: Real, M <: AbstractMeasure{X}}
    logdensity(sm.base, x) + sm.logscale
end

function Base.:*(k::Number, m::AbstractMeasure{X}) where {X} 
    if iszero(k)
        return TrivialMeasure(X)
    else
        return ScaledMeasure{typeof(k), typeof(m),X}(log(k),m)
    end
end

Base.:*(m::AbstractMeasure{X}, k::Real) where {X} = k * m

≪(::M, ::ScaledMeasure{R,M,X}) where {R,M,X} = true
≪(::ScaledMeasure{R,M,X}, ::M) where {R,M,X} = true
