struct ScaledMeasure{M,R}
    scale :: R
    base :: M
end

@implement IsMeasure{M,X} >: IsMeasure{ScaledMeasure{M,R},X} where {M, X, R <: Real}

@implement HasDensity{M,X} >: HasDensity{ScaledMeasure{M,R},X} where {M, X, R <: Real} begin
    logdensity(SM, X) = logdensity(SM.base, X) + log(SM.scale)
end



Base.:*(m::Measure{X}, k::Real) where {X} = ScaledMeasure{typeof(m),typeof(k)}(k,m)
Base.:*(k,m::Measure{X}) where {X} = m * k