# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: AbstractMeasure end

function Base.show(io::IO, μ::Lebesgue{X}) where X
    io = IOContext(io, :compact => true)
    print(io, "Lebesgue(", X, ")")
end

Lebesgue(X) = Lebesgue{X}()

basemeasure(μ::Lebesgue) = μ

isprimitive(::Lebesgue) = true

sampletype(::Lebesgue{ℝ}) = Float64

logdensity(::Lebesgue, x) = zero(float(x))
