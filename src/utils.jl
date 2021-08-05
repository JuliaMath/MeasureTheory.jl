const EmptyNamedTuple = NamedTuple{(),Tuple{}}

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

function fix(f, x)
    y = f(x)
    while x ≠ y
        (x,y) = (y, f(y))
    end

    return y
end

# function constructorof(::Type{T}) where {T} 
#     C = T
#     while C isa UnionAll
#         C = C.body
#     end

#     return C.name.wrapper
# end

constructor(::T) where {T} = constructor(T)

macro trysupport(ex)
    ex = esc(ex)
    quote
        try
            return $ex
        catch 
            return -Inf
        end
    end
end

export testvalue
testvalue(μ::AbstractMeasure) = testvalue(basemeasure(μ))

export rootmeasure

"""
    rootmeasure(μ::AbstractMeasure)

It's sometimes important to be able to find the fix point of a measure under
`basemeasure`. That is, to start with some measure and apply `basemeasure`
repeatedly until there's no change. That's what this does.
"""
rootmeasure(μ::AbstractMeasure) = fix(basemeasure, μ)

# Base on the Tricks.jl README
using Tricks
struct Iterable end
struct NonIterable end
isiterable(::Type{T}) where T = static_hasmethod(iterate, Tuple{T}) ? Iterable() : NonIterable()

functioninstance(::Type{F}) where {F<:Function} = F.instance

using Static
using ArrayInterface
using StaticArrays

struct KnownSize{S, T}
    value::T
end

KnownSize(x::T) where {T} = KnownSize{Tuple{ArrayInterface.known_size(T)...}, T}(x)
