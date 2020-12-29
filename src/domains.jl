abstract type AbstractDomain end

struct RealNumbers <: AbstractDomain end
ℝ = RealNumbers()
Base.show(io::IO, ::RealNumbers) = print(io, "ℝ")

struct PositiveReals <: AbstractDomain end
ℝ₊ = PositiveReals()
Base.show(io::IO, ::PositiveReals) = print(io, "ℝ₊")

struct UnitInterval <: AbstractDomain end
𝕀 = UnitInterval()
Base.show(io::IO, ::RealNumbers) = print(io, "𝕀")

struct Integers <: AbstractDomain end
ℤ = Integers()
Base.show(io::IO, ::Integers) = print(io, "ℤ")

struct NonnegativeIntegers <: AbstractDomain end
ℤ₊ = NonnegativeIntegers()
Base.show(io::IO, ::NonnegativeIntegers) = print(io, "ℤ₊")

struct PositiveIntegers <: AbstractDomain end
ℤ₀₊ = PositiveIntegers()
Base.show(io::IO, ::PositiveIntegers) = print(io, "ℤ₀₊")