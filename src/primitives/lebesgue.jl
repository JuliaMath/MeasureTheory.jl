# Lebesgue measure

export Lebesgue

struct Lebesgue{X} <: PrimitiveMeasure end

function Base.show(io::IO, μ::Lebesgue{X}) where X
    io = IOContext(io, :compact => true)
    print(io, "Lebesgue(")
    print(io, X)
    print(io, ")")
end

Lebesgue(X) = Lebesgue{X}()


sampletype(::Lebesgue{ℝ}) = Float64
sampletype(::Lebesgue{ℝ₊}) = Float64
sampletype(::Lebesgue{𝕀}) = Float64

testvalue(::Lebesgue{ℝ}) = 0.0
testvalue(::Lebesgue{𝕀}) = 0.5
testvalue(::Lebesgue{ℝ₊}) = 1.0
testvalue(::Lebesgue{<:Real}) = 0.0

logdensity(::Lebesgue, x) = zero(x)



Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = Lebesgue
