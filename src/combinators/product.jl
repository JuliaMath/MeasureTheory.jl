function TV.as(d::ProductMeasure{A}) where {A<:AbstractArray}
    d1 = marginals(d).f(first(marginals(d).data))
    as(Array, as(d1), size(marginals(d))...)
end

###############################################################################
# I <: Base.Generator

function TV.as(d::ProductMeasure{<:Base.Generator})
    d1 = marginals(d).f(first(marginals(d).iter))
    as(Array, as(d1), size(marginals(d))...)
end

function Base.rand(
    rng::AbstractRNG,
    ::Type,
    d::ProductMeasure{A},
    d1::Dists.Distribution,
) where {A<:AbstractArray}
    mar = marginals(d)

    # Distributions doens't (yet) have the three-argument form
    elT = typeof(rand(rng, first(mar)))

    sz = size(mar)
    x = Array{elT,length(sz)}(undef, sz)
    for (j, parj) in enumerate(mar)
        x[j] = rand(rng, d.f(parj))
    end
    x
end

# e.g. set(Normal(μ=2)^5, params, randn(5))
function Accessors.set(
    d::ProductMeasure{A},
    ::typeof(params),
    p::AbstractArray,
) where {A<:AbstractArray}
    set.(marginals(d), params, p)
end

function Accessors.set(
    d::ProductMeasure{A},
    ::typeof(params),
    p,
) where {M,A<:AbstractArray{M}}
    par = typeof(d.xs[1])(p)
    ProductMeasure(d.f, Fill(par, size(d.xs)))
end
