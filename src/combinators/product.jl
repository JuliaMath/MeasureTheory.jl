

using MappedArrays

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
