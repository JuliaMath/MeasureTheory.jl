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

function inferAdd(::Type{SuperpositionMeasure}, ::Type{SuperpositionMeasure})
    SuperpositionMeasure

function inferAdd end
@trait Add{A, B, C} where {C = inferAdd(A, B)} begin
    (+) :: [A, B] => C
    (+) = Base.:+
end

addtype()

@implement Add{}


function Base.:+(ms::(M where  {X, Measure{M, X} })...) 

    X = promote_type(eltype.(ms))
    measuretuples = getproperty.(ensureSuperposed.(ms), :components)
    SuperpositionMeasure(foldl(append!!, measuretuples))
end

