struct NegativeBinomial{N,T} <: ParameterizedMeasure{N}
    par :: NamedTuple{N,T}
end