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

Base.:*(k::Number, m::AbstractMeasure{X}) where {X} = ScaledMeasure{typeof(k), typeof(m),X}(log(k),m)

Base.:*(m::AbstractMeasure{X}, k::Real) where {X} = k * m
