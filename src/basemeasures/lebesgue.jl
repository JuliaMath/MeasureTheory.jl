# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: AbstractMeasure{X} end

Lebesgue(X) = Lebesgue{X}()

baseMeasure(μ::Lebesgue{X}) where {X} = μ
