export CountingMeasure

struct CountingMeasure{X} <: AbstractMeasure end

function Base.show(io::IO, μ::CountingMeasure{X}) where {X}
    io = IOContext(io, :compact => true)
    print(io, "CountingMeasure(", X, ")")
end

CountingMeasure(X) = CountingMeasure{X}()

basemeasure(μ::CountingMeasure) = μ

isprimitive(::CountingMeasure) = true

# sampletype(::CountingMeasure{ℝ}) = Float64
# sampletype(::CountingMeasure{ℝ₊}) = Float64
# sampletype(::CountingMeasure{𝕀}) = Float64

sampletype(::CountingMeasure{IntegerRange{lo,hi}}) where {lo, hi} = Int



logdensity(::CountingMeasure, x) = zero(float(x))

# (::CountingMeaure)(s) = length(Set(s))
