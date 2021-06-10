

using SimpleTraits


# Type hierarchy for cheap implications
abstract type MeasureTrait{X} <: SimpleTraits.Trait end
abstract type IsRepresentative{X} <: MeasureTrait{X} end
abstract type IsPrimitive{X} <: IsRepresentative{X} end

abstract type IsScalar{X} <: MeasureTrait{X} end

"""
    @scalar T

Declare that every measure of type `T` are scalar.
"""
macro scalar(T)
    MT = MeasureTheory
    ST = MeasureTheory.SimpleTraits

    # Really the same as `@traitimpl IsScalar{T}`
    # But maybe easier for beginners
    quote
        function $ST.trait(::$ST.Type{$MT.IsScalar{S}}) where {S <: $T}
            return $MT.IsScalar{S} 
        end
    end
end

"""
    isscalar(::AbstractMeasure)
    isscalar(::Type{M}) where {M<:AbstractMeasure} 

Tests whether a given measure has been declared scalar. See `@scalar`.
"""
function isscalar end

isscalar(::M) where {M} = isscalar(M)
isscalar(::Type{M}) where {M} = istrait(IsScalar{T})

"""
    isprimitive(::AbstractMeasure)
    isprimitive(::Type{M}) where {M<:AbstractMeasure} 

Most measures are defined in terms of other measures, for example using a
density or a pushforward. Those that are not are considered (in this library,
it's not a general measure theory thing) to be _primitive_. The canonical
example of a primitive measure is `Lebesgue(X)` for some `X`.
"""
function isprimitive end

isprimitive(::M) where {M} = isprimitive(M)
isprimitive(::Type{M}) where {M} = istrait(IsPrimitive{T})

"""
    @primitive T

Declare that every measure of type `T` should be considered "primitive". A
measure is primitive is it is not defined in terms of another measure. Common
examples are Lebesgue and Counting measures. 

Note that this is not a general measure-theoretic term, but is specific to the
MeasureTheory.jl implementation. 
"""
macro primitive(T)
    MT = MeasureTheory
    ST = MeasureTheory.SimpleTraits

    # Really the same as `@traitimpl IsPrimitive{T}`
    # But maybe easier for beginners
    quote
        function $ST.trait(::$ST.Type{$MT.IsPrimitive{S}}) where {S <: $T}
            return $MT.IsPrimitive{S} 
        end
    end
end


export basemeasure

"""
    basemeasure(μ)

Many measures are defined in terms of a logdensity relative to some base
measure. This makes it important to be able to find that base measure.

For measures not defined in this way, we'll typically have `basemeasure(μ) == μ`.
"""
function basemeasure end
