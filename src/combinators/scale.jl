struct ScaledMeasure{R,M,X} <: Measure{X}
    scale :: R
    base :: M
end

function logdensity(sm::ScaledMeasure{R,M,X}, x::X) where {X, R <: Real, M <: Measure{X}}
    logdensity(sm.base, x) + log(sm.scale)
end

Base.:*(m::Measure{X}, k::Real) where {X} = ScaledMeasure{typeof(k), typeof(m),X}(k,m)
Base.:*(k::Real, m::Measure{X}) where {X} = m * k
