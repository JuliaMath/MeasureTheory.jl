struct ScaledMeasure{X,R} <: Measure{X}
    scale :: R
    base :: Measure{X}
end

Base.:*(m::Measure{X}, k::Real) where {X} = ScaledMeasure{X,typeof(k)}(k,m)
Base.:*(k,m::Measure{X}) where {X} = m * k