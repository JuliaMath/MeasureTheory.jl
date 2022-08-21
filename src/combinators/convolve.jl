struct Convolution{M,N} <: AbstractMeasure
    μ::M
    ν::N
end

"""
If μ, ν are subtypes of `AbstractMeasure` or satisfy the Measure interface, r
then `convolve(μ, ν)` is a measure, called the convolution of μ and ν.
"""
convolve(μ, ν) = Convolution(μ, ν)

function Base.rand(rng::AbstractRNG, ::Type{T}, d::Convolution) where {T}
    x = rand(rng, T, d.μ)
    y = rand(rng, T, d.ν)
    return x+y
end
