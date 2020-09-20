export SuperpositionMeasure

struct SuperpositionMeasure{X,NT} <: AbstractMeasure{X}
    components :: NT
end

SuperpositionMeasure(ms :: AbstractMeasure{X}...) where {X} = SuperpositionMeasure{X,length(ms)}(ms)

# SuperpositionMeasure(m::NTuple{N, Measure{X}}) where {N,X} = SuperpositionMeasure(m...)

# Base.length(::SuperpositionMeasure{X,N}) where {X,N} = N

function Base.:+(μ::SuperpositionMeasure{X,N1}, ν::SuperpositionMeasure{X,N2}) where {X,N1,N2}
    components = (μ.components..., ν.components...)
    SuperpositionMeasure{X, N1+N2}(components)
end

function Base.:+(μ::AbstractMeasure{X}, ν::SuperpositionMeasure{X,N}) where {X,N}
    components = (μ, ν.components...)
    SuperpositionMeasure{X,N+1}(components)
end

function Base.:+(μ::SuperpositionMeasure{X,N}, ν::AbstractMeasure{X}) where {X,N}
    components = (μ.components..., ν)
    SuperpositionMeasure{X,N+1}(components)
end

function Base.:+(μ::AbstractMeasure{X}, ν::AbstractMeasure{X}) where {X}
    components = (μ, ν)
    SuperpositionMeasure{X,2}(components)
end

# TODO : rand(::SuperpositionMeasure)
# function Base.rand(μ::SuperpositionMeasure{X,N}) where {X,N}
# end
