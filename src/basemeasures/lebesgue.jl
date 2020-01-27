# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: Measure{X} end

@implement IsMeasure{Lebesgue{X}, X} where {X <: Real}

@implement IsMeasure{Lebesgue{X}, X} where {X <: Vector{Real}}

Lebesgue(X) = Lebesgue{X}()