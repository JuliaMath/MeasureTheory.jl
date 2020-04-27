export SuperpositionMeasure

struct SuperpositionMeasure{X,N} <: Measure{X}
    components :: NTuple{N,Measure{X}}
end

SuperpositionMeasure(ms :: Measure{X}...) where {X} = SuperpositionMeasure{X,length(ms)}(ms)

# SuperpositionMeasure(m::NTuple{N, Measure{X}}) where {N,X} = SuperpositionMeasure(m...)

Base.length(::SuperpositionMeasure{X,N}) where {X,N} = N

function +(μ::SuperpositionMeasure{X,N1}, ν::SuperpositionMeasure{X,N2}) where {X,N1,N2}
    components = append!!(μ.components, ν.components)
    SuperpositionMeasure{X, N1+N2}(components)
end
