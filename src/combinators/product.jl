function as(d::PowerMeasure)
    as(Array, as(d.parent), length.(d.axes)...)
end

function as(d::ProductMeasure{<:AbstractArray{<:Dirac}})
    return asConst(testvalue.(marginals(d)))
end

function as(d::ProductMeasure{A}) where {A<:AbstractArray}
    mar = marginals(d)
    ts = map(as, mar)
    if allequal(ts)
        return as(Array, first(ts), size(ts))
    else
        error("Not yet implemented")
    end
end

using MappedArrays

function as(d::ProductMeasure{A}) where {A<:MappedArrays.ReadonlyMappedArray}
    d1 = marginals(d).f(first(marginals(d).data))
    as(Array, as(d1), size(marginals(d))...)
end

function as(d::ProductMeasure{T}) where {T<:Tuple}
    as(map(as, d.marginals))
end

###############################################################################
# I <: Base.Generator

function as(d::ProductMeasure{<:Base.Generator})
    d1 = marginals(d).f(first(marginals(d).iter))
    as(Array, as(d1), size(marginals(d))...)
end

# function as(d::ProductMeasure{Returns{T},F,A}) where {T,F,A<:AbstractArray}
#     as(Array, as(d.f.f.value), size(d.xs))
# end

@inline function Base.rand(
    rng::AbstractRNG,
    ::Type{T},
    d::ProductMeasure{A},
) where {T,A<:AbstractArray}
    mar = marginals(d)

    # Distributions doens't (yet) have the three-argument form
    elT = typeof(rand(rng, T, first(mar)))

    sz = size(mar)
    x = Array{elT,length(sz)}(undef, sz)
    @inbounds @simd for j in eachindex(mar)
        x[j] = rand(rng, T, mar[j])
    end
    x
end

# # e.g. set(Normal(μ=2)^5, params, randn(5))
# function Accessors.set(
#     d::ProductMeasure{A},
#     ::typeof(params),
#     p::AbstractArray,
# ) where {A<:AbstractArray}
#     set.(marginals(d), params, p)
# end

# function Accessors.set(
#     d::ProductMeasure{A},
#     ::typeof(params),
#     p,
# ) where {A<:AbstractArray}
#     mar = marginals(d)
#     par = eltype(mar)(p)
#     ProductMeasure(d.f, Fill(par, size(mar)))
# end
