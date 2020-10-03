# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: AbstractMeasure end

Lebesgue(X) = Lebesgue{X}()

basemeasure(μ::Lebesgue{X}) where {X} = μ

isprimitive(::Lebesgue) = true

sampletype(::Lebesgue{X}) where{X} = X
