
export For
using Random
import Base

For{D,N,T,F} = ProductMeasure{ReadonlyMappedArray{D, N, T, F}} 

function Base.show(io::IO, d::ProductMeasure{<:AbstractMappedArray})
    print(io, "For(")
    print(io, d.data.f, ", ")
    print(io, d.data.data, ")")
end

For(f, dims) = ProductMeasure(mappedarray(f, dims))


function Base.eltype(::For{D,N,T,F}) where {D,N,T,F}
    return eltype(D)
end


# """
#     indexstyle(a::AbstractArray, b::AbstractArray)

# Find the best IndexStyle that works for both `a` and `b`. This will return
# `IndexLinear` if both `a` and `b` support it; otherwise it will fall back on `IndexCartesian`.
# """
# function indexstyle(::A,::B)
#     if IndexStyle(A) == IndexStyle(B) == IndexLinear()
#         return IndexLinear()
#     end

#     return IndexCartesian()
# end

# function Base.rand(rng::AbstractRNG, μ::For{D,N,T,F}) where {F,T<:AbstractArray,D,X}
#     s = size(μ.θ)
#     x = Array{X,length(s)}(undef, s...)
#     rand!(rng, x, μ)
# end

# function logdensity(μ::For{D,N,T,F}, x)
#     getℓ(θⱼ, xⱼ) = logdensity(μ.f(θⱼ), xⱼ)
#     ℓ = mappedarray(getℓ, μ.θ, x)
#     _logdensity(μ, x, indexstyle(μ.θ, x), result_type)
# end

# function _logdensity(μ::For{D,N,T,F}, x, ::IndexLinear, ::Type{R}) where {R<:AbstractFloat}
#     ℓ = zero(R)
#     μ.f(μ.θ)
# end

# function basemeasure(μ::For{D,N,T,F}) where {F,T<:AbstractArray,D,X}
