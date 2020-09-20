# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: AbstractMeasure{X} end

Lebesgue(X::DataType) = Lebesgue{X}()

add!(BASE_MEASURES, Lebesgue{X} where X)
