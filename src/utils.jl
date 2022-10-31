using LinearAlgebra

function solve(A::Union{AbstractMatrix,Factorization}, y::AbstractArray)
    (m, n) = size(A)
    n == 1 && return [dot(A, y) / sum(a -> a^2, A)]
    return A \ y
end

@inline function mydot(a::Real, b::Real)
    return a * b
end

@inline function mydot(a::Ta, b::Tb) where {Ta<:AbstractArray,Tb<:AbstractArray}
    return dot(a, b)
end

@inline function _mydot(a::NTuple{0}, b::NTuple{0}, s)
    return s
end

@inline function mydot(a::Tuple, b::Tuple)
    z = 0.0
    _mydot(a, b, z)
end

@inline function _mydot(a::Tuple, b::Tuple, s)
    _mydot(Base.tail(a), Base.tail(b), s + mydot(first(a), first(b)))
end

@generated function mydot(a::SVector{N,Ta}, b::SVector{N,Tb}) where {N,Ta,Tb}
    z = zero(float(Base.promote_eltype(Ta, Tb)))
    quote
        $(Expr(:meta, :inline))
        result = $z
        @inbounds Base.Cartesian.@nexprs $N i -> begin
            result += a[i] * b[i]
        end
        return result
    end
end

in𝕀(x) = static(0.0) ≤ x ≤ static(1.0)
inℝ₊(x) = static(0.0) ≤ x

_eltype(T) = eltype(T)

function _eltype(g::Base.Generator{I}) where {I}
    Core.Compiler.return_type(g.f, Tuple{_eltype(I)})
end

function _eltype(
    ::Type{Base.Generator{I,ComposedFunction{Outer,Inner}}},
) where {Outer,Inner,I}
    _eltype(Base.Generator{_eltype(Base.Generator{I,Inner}),Outer})
end

function _eltype(::Type{Base.Generator{I,F}}) where {F<:Function,I}
    f = instance(F)
    Core.Compiler.return_type(f, Tuple{_eltype(I)})
end

function _eltype(::Type{Z}) where {Z<:Iterators.Zip}
    map(_eltype, Z.types[1].types)
end

# Adapted from https://github.com/JuliaArrays/MappedArrays.jl/blob/46bf47f3388d011419fe43404214c1b7a44a49cc/src/MappedArrays.jl#L229
function func_string(f, types)
    ft = typeof(f)
    mt = ft.name.mt
    name = string(mt.name)
    if startswith(name, '#')
        # This is an anonymous function. See if it can be printed nicely
        lwrds = code_lowered(f, types)
        if length(lwrds) != 1
            return string(f)
        end
        lwrd = lwrds[1]
        c = lwrd.code
        if length(c) == 2 && (
            (isa(c[2], Expr) && c[2].head === :return) ||
            (isdefined(Core, :ReturnNode) && isa(c[2], Core.ReturnNode))
        )
            # This is a single-line anonymous function, we should handle it
            s = lwrd.slotnames[2:end]
            result = ""
            if length(s) == 1
                result *= string(s[1])
                result *= "->"
            else
                result *= "(" * join(string.(tuple(s...)), ", ") * ")"
                result *= "->"
            end
            c1 = string(c[1])
            for i in 1:length(s)
                c1 = replace(c1, "_" * string(i + 1) => string(s[i]))
            end
            result *= c1
            return result
        else
            return string(f)
        end
    else
        return string(f)
    end
end

function getL(C::Cholesky)
    Cfactors = getfield(C, :factors)
    Cuplo    = getfield(C, :uplo)

    LowerTriangular(Cuplo === 'L' ? Cfactors : Cfactors')
end

function getU(C::Cholesky)
    Cfactors = getfield(C, :factors)
    Cuplo    = getfield(C, :uplo)

    UpperTriangular(Cuplo === 'U' ? Cfactors : Cfactors')
end

const Triangular = Union{L,U} where {L<:LowerTriangular,U<:UpperTriangular}

# Ideally we could just add one method smf(μ::$M, x::Dual{TAG}) where TAG
# But then we have method ambiguities
macro smfAD(M)
    quote
        function MeasureBase.smf(μ::$M, x::Dual{TAG}) where {TAG}
            val = ForwardDiff.value(x)
            Δ = ForwardDiff.partials(x)
            Dual{TAG}(smf(μ, val), Δ * densityof(μ, val))
        end
    end
end
