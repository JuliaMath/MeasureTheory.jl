export ProductMeasure

struct ProductMeasure{T} <: Measure{T}
    # e.g. Tuple{Int64,Float64} => Tuple{Measure{Int64},Measure{Float64}}
    components :: Tuple{(Measure{t} for t in T.types)...}
end

ProductMeasure

# ProductMeasure(m::NTuple{N, Measure{X}}) where {N,X} = ProductMeasure(m...)

Base.length(m::ProductMeasure{T}) where {T} = length(m)


function Base.:*(μ::ProductMeasure{X,N1}, ν::ProductMeasure{X,N2})
    components = append!!(μ.components, ν.components)
    ProductMeasure{X, N1+N2}(components)
end
