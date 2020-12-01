
export For
using Random

struct For{F,T,D,X}
    f::F
    θ::T
end

function Base.eltype(::For{F,T,D,X}) where {F,T,D,X}
    return X
end


"""
    indexstyle(a::AbstractArray, b::AbstractArray)

Find the best IndexStyle that works for both `a` and `b`. This will return
`IndexLinear` if both `a` and `b` support it; otherwise it will fall back on `IndexCartesian`.
"""
@generated function indexstyle(a,b)
    if IndexStyle(a) == IndexStyle(b) == IndexLinear()
        return IndexLinear()
    end

    return IndexCartesian()
end

Base.rand(μ::For) = rand(GLOBAL_RNG, dist)

import Random

function Random.rand!(rng::AbstractRNG, x::AbstractArray, μ::For{F,T,D,X}) where {F,T<:AbstractArray,D,X}
    @inbounds @simd for i in eachindex(indexstyle(x, μ.θ), x)
        # TODO: What if `μ.f(μ.θ[i])` allocates? Can we use `rand!` here too?
        x[i] = rand(rng, μ.f(μ.θ[i]))
    end
end

# function basemeasure(μ::For{F,T,D,X}) where {F,T<:AbstractArray,D,X}
