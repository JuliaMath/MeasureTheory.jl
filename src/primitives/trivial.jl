export TrivialMeasure

struct TrivialMeasure <: PrimitiveMeasure end

@primitive TrivialMeasure

sampletype(::TrivialMeasure) = Nothing
