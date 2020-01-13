export SuperpositionMeasure

struct SuperpositionMeasure{X,N} <: Measure{X}
    components :: NTuple{N,Measure{X}}
    SuperpositionMeasure(ms :: Measure{X} ...) where {X} = new{X,length(ms)}(ms)
end

SuperpositionMeasure(m::NTuple{N, Measure{X}}) where {N,X} = SuperpositionMeasure(m...)

Base.length(::SuperpositionMeasure{X,N}) where {X,N} = N

ensureSuperposed(m::SuperpositionMeasure{X,N}) where {X,N} = m
ensureSuperposed(m::Measure{X}) where {X} = SuperpositionMeasure(m)

using BangBang

function Base.:+(ms::Measure{X}...) where {X}
    measuretuples = getproperty.(ensureSuperposed.(ms), :components)
    SuperpositionMeasure(foldl(append!!, measuretuples))
end
