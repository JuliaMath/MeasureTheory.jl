# abstract type AbstractOrderedVector{T} <: AbstractVector{T} end

# struct OrderedVector{T} <: AbstractOrderedVector{T}
#     data::Vector{T}
# end

# as(::Type{OrderedVector}, transformation::TV.AbstractTransform, dim::Int) =
#     Ordered(transformation, dim)

using MeasureTheory

export Ordered

using TransformVariables
const TV = TransformVariables

struct Ordered{T <: TV.AbstractTransform} <: TV.VectorTransform
    transformation::T
    dim::Int
end

TV.dimension(t::Ordered) = t.dim

addlogjac(ℓ, Δℓ) = ℓ + Δℓ
addlogjac(::TV.NoLogJac, Δℓ) = TV.NoLogJac()

using MappedArrays

bounds(t::TV.ShiftedExp{true}) = (t.shift, TV.∞)
bounds(t::TV.ShiftedExp{false}) = (-TV.∞, t.shift)
bounds(t::TV.ScaledShiftedLogistic) = (t.shift, t.scale + t.shift)
bounds(::TV.Identity) = (-TV.∞, TV.∞)

const OrderedΔx = -8.0

# See https://mc-stan.org/docs/2_27/reference-manual/ordered-vector.html
function TV.transform_with(flag::TV.LogJacFlag, t::Ordered, x, index::T) where {T}
    transformation, len = (t.transformation, t.dim)
    @assert dimension(transformation) == 1
    y = similar(x, len)
        
    (lo,hi) = bounds(transformation)


    x = mappedarray(xj -> xj + OrderedΔx, x)

    @inbounds (y[1], ℓ, _) = TV.transform_with(flag, as(Real, lo, hi), x, index)
    index += 1

    @inbounds for i in 2:len
        (y[i], Δℓ, _) =  TV.transform_with(flag, as(Real, y[i-1], hi), x, index)
        ℓ = addlogjac(ℓ, Δℓ)
        index += 1
    end

    return (y, ℓ, index)
end

TV.inverse_eltype(t::Ordered, y::AbstractVector) = TV.extended_eltype(y)

Ordered(n::Int) = Ordered(asℝ, n)

function TV.inverse_at!(x::AbstractVector, index, t::Ordered, y::AbstractVector)
    (lo,hi) = bounds(t.transformation)

    @inbounds x[index] = inverse(as(Real, lo, hi), y[1]) - OrderedΔx
    index += 1

    @inbounds for i in 2:length(y)
        x[index] = inverse(as(Real, y[i-1], hi), y[i]) - OrderedΔx
        index += 1
    end
    return index
end

export Sorted
struct Sorted{M} <: AbstractMeasure
    μ::M
    n::Int
end

logdensity(s::Sorted, x) = logdensity(s.μ ^ s.n, x)

TV.as(s::Sorted) = Ordered(as(s.μ), s.n)

function Random.rand!(rng::AbstractRNG, d::Sorted, x::AbstractArray)
    rand!(rng, d.μ ^ d.n, x)
    sort!(x)
    return x
end

function Base.rand(rng::AbstractRNG, T::Type, d::Sorted)
    # TODO: Use `sampletype` for this
    elT = typeof(rand(rng, T, d.μ))
    x = Vector{elT}(undef, d.n)
    rand!(rng, d, x)
end

# logdensity(d, rand(d))



# TV.transform_with(TV.LogJac(), Ordered(asℝ, 4), zeros(4), 1)
# TV.transform_with(TV.LogJac(), Ordered(asℝ, 4), randn(4), 1)

# d = Pushforward(Ordered(10), Normal()^10, false)
# logdensity(Lebesgue(ℝ)^10, rand(d))
