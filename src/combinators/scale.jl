struct ScaledMeasure{R,M,X} <: AbstractMeasure{X}
    scale :: R
    base :: M
end

function logdensity(sm::ScaledMeasure{R,M,X}, x::X) where {X, R <: Real, M <: AbstractMeasure{X}}
    logdensity(sm.base, x) + log(sm.scale)
end

Base.:*(m::AbstractMeasure{X}, k::Real) where {X} = ScaledMeasure{typeof(k), typeof(m),X}(k,m)
Base.:*(k::Real, m::AbstractMeasure{X}) where {X} = m * k
