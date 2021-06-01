export TrivialMeasure

struct TrivialMeasure <: AbstractMeasure end

@primitive TrivialMeasure

sampletype(::TrivialMeasure) = Nothing
