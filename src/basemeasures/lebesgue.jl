# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: AbstractMeasure end

function Base.show(io::IO, μ::Lebesgue{X}) where X
    io = IOContext(io, :compact => true)
    print(io, "Lebesgue(", X, ")")
end

Lebesgue(X) = Lebesgue{X}()

basemeasure(μ::Lebesgue{X}) where {X} = μ

isprimitive(::Lebesgue) = true

sampletype(::Lebesgue{X}) where{X} = X
