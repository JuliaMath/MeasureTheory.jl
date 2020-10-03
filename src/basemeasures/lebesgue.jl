# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: AbstractMeasure end

Lebesgue(X) = Lebesgue{X}()

baseMeasure(μ::Lebesgue{X}) where {X} = μ

isprimitive(::Lebesgue) = true
isprimitive(μ) = false

sampletype(::Lebesgue{X}) where{X} = X
