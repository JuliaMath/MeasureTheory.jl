# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: Measure{X} end
Lebesgue(X) = Lebesgue{X}()