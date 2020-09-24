struct TrivialMeasure{X} <: AbstractMeasure{X} end

TrivialMeasure(X) = TrivialMeasure{X}()
