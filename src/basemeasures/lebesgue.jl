# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: AbstractMeasure{X} end

Lebesgue(X::DataType) = Lebesgue{X}()
