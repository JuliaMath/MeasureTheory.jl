export SuperpositionMeasure

struct SuperpositionMeasure{X,N} <: Measure{X}
    components :: NTuple{N,Measure{X}}
end

SuperpositionMeasure(ms :: Measure{X}...) where {X} = SuperpositionMeasure{X,length(ms)}(ms)

# SuperpositionMeasure(m::NTuple{N, Measure{X}}) where {N,X} = SuperpositionMeasure(m...)

# Base.length(::SuperpositionMeasure{X,N}) where {X,N} = N

function Base.:+(μ::SuperpositionMeasure{X,N1}, ν::SuperpositionMeasure{X,N2}) where {X,N1,N2}
    components = (μ.components..., ν.components...)
    SuperpositionMeasure{X, N1+N2}(components)
end

function Base.:+(μ::Measure{X}, ν::SuperpositionMeasure{X,N}) where {X,N}
    components = (μ, ν.components...)
    SuperpositionMeasure{X,N+1}(components)
end

function Base.:+(μ::SuperpositionMeasure{X,N}, ν::Measure{X}) where {X,N}
    components = (μ.components..., ν)
    SuperpositionMeasure{X,N+1}(components)
end

function Base.:+(μ::Measure{X}, ν::Measure{X}) where {X}
    components = (μ, ν)
    SuperpositionMeasure{X,2}(components)
end

# TODO : rand(::SuperpositionMeasure)
# function Base.rand(μ::SuperpositionMeasure{X,N}) where {X,N}
# end
