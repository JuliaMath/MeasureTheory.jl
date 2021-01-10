abstract type AbstractDomain end

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
