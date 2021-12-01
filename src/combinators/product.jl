function TV.as(d::ProductMeasure{F,S,A}) where {F,S,A<:AbstractArray}
    d1 = marginals(d).f(first(marginals(d).data))
    as(Array, as(d1), size(marginals(d))...)
end

###############################################################################
# I <: Base.Generator

function TV.as(d::ProductMeasure{K,I}) where {K,I<:Base.Generator}
    d1 = marginals(d).f(first(marginals(d).iter))
    as(Array, as(d1), size(marginals(d))...)
end

function TV.as(d::ProductMeasure{Returns{T},F,A}) where {T,F,A<:AbstractArray}
    as(Array, as(d.f.f.value), size(d.xs))
end

function Base.rand(
    rng::AbstractRNG,
    ::Type{T},
    d::ProductMeasure,
    d1::Dists.Distribution,
) where {T}
    mar = marginals(d)

    # Distributions doens't (yet) have the three-argument form
    elT = typeof(rand(rng, first(mar)))

    sz = size(mar)
    x = Array{elT,length(sz)}(undef, sz)
    for (j, parj) in enumerate(d.xs)
        x[j] = rand(rng, d.f(parj))
    end
    x
end

# e.g. set(Normal(μ=2)^5, params, randn(5))
function Accessors.set(
    d::ProductMeasure{F,S,A},
    ::typeof(params),
    p::AbstractArray,
) where {F,S,A<:AbstractArray}
    set.(marginals(d), params, p)
end

function Accessors.set(
    d::ProductMeasure{F,S,A},
    ::typeof(params),
    p,
) where {F,S,A<:AbstractArray}
    par = typeof(d.xs[1])(p)
    ProductMeasure(d.f, Fill(par, size(d.xs)))
end
