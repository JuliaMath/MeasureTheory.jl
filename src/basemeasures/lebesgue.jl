# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: Measure{X} end

Lebesgue(X::DataType) = Lebesgue{X}()

add!(BASE_MEASURES, Lebesgue{X} where X)
