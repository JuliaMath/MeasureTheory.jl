export SuperpositionMeasure

struct SuperpositionMeasure{X,N} <: Measure{X}
    components :: NTuple{N,Measure{X}}
end

SuperpositionMeasure(ms :: Measure{X}...) where {X} = SuperpositionMeasure{X,length(ms)}(ms)

# SuperpositionMeasure(m::NTuple{N, Measure{X}}) where {N,X} = SuperpositionMeasure(m...)

Base.length(::SuperpositionMeasure{X,N}) where {X,N} = N

ensureSuperposed(m::SuperpositionMeasure{X,N}) where {X,N} = m
ensureSuperposed(m::Measure{X}) where {X} = SuperpositionMeasure(m)

using BangBang

function inferAdd(::Type{SuperpositionMeasure{X,N1}}, ::Type{SuperpositionMeasure{X,N2}}) where {X,N1,N2}
    SuperpositionMeasure{X, N1+N2}
end

function inferAdd end
@trait Add{A, B} begin
    (+) :: [A, B] => inferAdd(A, B)
    (+) = Base.:+
end

export Add

@implement Add{SuperpositionMeasure{X,N1}, SuperpositionMeasure{X,N2}} where {X, N1, N2} begin
    function +(μ,ν) 
        components = append!!(μ.components, ν.components)
        SuperpositionMeasure{X, N1+N2}(components)
    end
end


