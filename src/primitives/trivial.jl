export TrivialMeasure

struct TrivialMeasure <: PrimitiveMeasure end

sampletype(::TrivialMeasure) = Nothing
