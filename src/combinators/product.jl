export ProductMeasure

struct ProductMeasure{T} <: Measure{T}
    # e.g. Tuple{Int64,Float64} => Tuple{Measure{Int64},Measure{Float64}}
    components :: Tuple{(Measure{t} for t in T.types)...}
end

ProductMeasure

# ProductMeasure(m::NTuple{N, Measure{X}}) where {N,X} = ProductMeasure(m...)

Base.length(m::ProductMeasure{T}) where {T} = length(m)


function inferMul(::Type{ProductMeasure{X,N1}}, ::Type{ProductMeasure{X,N2}}) where {X,N1,N2}
    ProductMeasure{X, N1+N2}
end

function inferMul end
@trait Mul{A, B} begin
    (*) :: [A, B] => inferMul(A, B)
    (*) = Base.:*
end

export Mul

@implement Mul{ProductMeasure{X,N1}, ProductMeasure{X,N2}} where {X, N1, N2} begin
    function *(μ,ν) 
        components = append!!(μ.components, ν.components)
        ProductMeasure{X, N1+N2}(components)
    end
end
