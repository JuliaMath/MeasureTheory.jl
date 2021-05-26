abstract type AbstractDomain end

"""
    @domain(name, T)

Defines a new singleton struct `T`, and a value `name` for building values of
that type.

For example, `@domain ℝ RealNumbers` is equivalent to

    struct RealNumbers <: AbstractDomain end

    export ℝ

    ℝ = RealNumbers()

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


###########################################################
# Integer ranges

struct IntegerRange{lo, hi} <: AbstractDomain end

Base.minimum(::IntegerRange{lo, hi}) where {lo, hi} = lo
Base.maximum(::IntegerRange{lo, hi}) where {lo, hi} = hi

iterate(r::IntegerRange{lo, hi}) where {lo, hi} = iterate(lo:hi)

function Base.getindex(::Integers, r::AbstractUnitRange)
    IntegerRange{minimum(r), maximum(r)}()
end

function Base.show(io::IO, r::IntegerRange{lo, hi}) where {lo, hi}
    print(io, "ℤ[", lo, ":", hi, "]")
end

testvalue(::IntegerRange{lo, hi}) where {lo, hi} = lo

###########################################################
# Real intervals


struct RealInterval{lo, hi} <: AbstractDomain end
