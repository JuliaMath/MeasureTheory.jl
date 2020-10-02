export TrivialMeasure

struct TrivialMeasure <: AbstractMeasure end

sampletype(::TrivialMeasure) = Nothing
