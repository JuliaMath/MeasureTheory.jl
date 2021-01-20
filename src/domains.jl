abstract type AbstractDomain end

"""
    @domain(name, T)

Defines a new singleton struct `T`, and a value `name` for building values of
that type.

For example, `MeasureTheory.@domain ℝ RealNumbers` is equivalent to

    struct RealNumbers <: MeasureTheory.AbstractDomain end

    export ℝ

    ℝ = MeasureTheory.RealNumbers()

    Base.show(io::IO, ::RealNumbers) = print(io, "ℝ")
"""
macro domain(name, T)
    sname = String(name)
    
    name = esc(name)
    quote
        struct $T <: AbstractDomain end
        export $name
        $name = $T()
        Base.show(io::IO, ::$T) = print(io, $sname)
    end
end

@domain ℝ RealNumbers

@domain ℝ₊ PositiveReals

@domain 𝕀 UnitInterval

@domain ℤ Integers

@domain ℤ₊ PositiveIntegers

@domain ℤ₀₊ NonnegativeIntegers
