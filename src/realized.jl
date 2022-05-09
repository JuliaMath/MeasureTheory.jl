using DynamicIterators
import DynamicIterators: dub, dyniterate, evolve

struct SamplesTo{M,T}
    measure::M
    element::T
end

export ↝

↝(m, x) = SamplesTo(m, x)

function Pretty.tile(s::SamplesTo)
    Pretty.pair_layout(Pretty.tile(s.measure), Pretty.tile(s.element); sep = " ↝ ")
end

function Base.show(io::IO, s::SamplesTo)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, s)
end

###############################################################################

struct Realized{R,S,T} <: DynamicIterators.DynamicIterator
    rng::ResettableRNG{R,S}
    iter::T
end

Base.copy(r::Realized) = Realized(copy(r.rng), r.iter)

reset!(rv::Realized) = reset!(rv.rng)

Base.show(io::IO, r::Realized) = Pretty.pprint(io, r)

function Pretty.quoteof(r::Realized)
    :(Realized($(Pretty.quoteof(r.rng)), $(Pretty.quoteof(r.iter))))
end

Base.IteratorEltype(mc::Realized) = Base.HasEltype()

# function Base.eltype(::Type{Rz}) where {R,S,T,Rz <: Realized{R,S,T}}
#     eltype(T)
# end

Base.length(r::Realized) = length(r.iter)

Base.size(r::Realized) = size(r.iter)

Base.IteratorSize(::Type{Rz}) where {R,S,T,Rz<:Realized{R,S,T}} = Base.IteratorSize(T)
Base.IteratorSize(r::Rz) where {R,S,T,Rz<:Realized{R,S,T}} = Base.IteratorSize(r.iter)

Base.iterate(rv::Realized) = iterate(rv, nothing)

function Base.iterate(rv::Realized{R,S,T}, ::Nothing) where {R,S,T}
    reset!(rv.rng)
    if static_hasmethod(evolve, Tuple{T})
        dyniterate(rv, nothing)
    else
        μs = iterate(rv.iter)
        isnothing(μs) && return nothing
        (μ, s) = μs
        x = rand(rv.rng, μ)
        return (μ ↝ x), s
    end
end

function Base.iterate(rv::Realized{R,S,T}, s) where {R,S,T}
    if static_hasmethod(evolve, Tuple{T})
        dyniterate(rv, s)
    else
        μs = iterate(rv.iter, s)
        isnothing(μs) && return nothing
        (μ, s) = μs
        x = rand(rv.rng, μ)
        return (μ ↝ x), s
    end
end

function dyniterate(rv::Realized, ::Nothing)
    rv = copy(rv)
    reset!(rv.rng)
    μ = evolve(rv.iter)
    x = rand(rv.rng, μ)
    (μ ↝ x), Sample(x, rv.rng)
end

function dyniterate(rv::Realized, u::Sample)
    dyniterate(rv.iter, u)
end

############################################################

struct RealizedMeasures{R,S,T} <: DynamicIterators.DynamicIterator
    rng::ResettableRNG{R,S}
    iter::T
end

Base.show(io::IO, r::RealizedMeasures) = Pretty.pprint(io, r)

function Pretty.quoteof(r::RealizedMeasures)
    :(RealizedMeasures($(Pretty.quoteof(r.rng)), $(Pretty.quoteof(r.iter))))
end

function Base.iterate(rm::RealizedMeasures, s = nothing)
    val, s = iterate(rm, s)
    (val.measure, s)
end

function dyniterate(rm::RealizedMeasures, s)
    val, s = dyniterate(rm, s)
    (val.measure, s)
end

##########################################

struct RealizedSamples{R,S,T} <: DynamicIterators.DynamicIterator
    parent::Realized{R,S,T}
end

RealizedSamples(rng::AbstractRNG, iter) = RealizedSamples(Realized(rng, iter))

Base.show(io::IO, r::RealizedSamples) = Pretty.pprint(io, r)

function Pretty.quoteof(r::RealizedSamples)
    :(RealizedSamples($(Pretty.quoteof(r.parent.rng)), $(Pretty.quoteof(r.parent.iter))))
end

function iter_helper(x)
    isnothing(x) && return x
    val, s = x
    (val.element, s)
end

function Base.iterate(rm::RealizedSamples)
    iter_helper(iterate(rm.parent))
end

function Base.iterate(rm::RealizedSamples, ::Nothing)
    iter_helper(iterate(rm.parent))
end

function Base.iterate(rm::RealizedSamples, s)
    iter_helper(iterate(rm.parent, s))
end

function dyniterate(rm::RealizedSamples, s)
    val, s = dyniterate(rm.parent, s)
    (val.element, s)
end

Base.copy(r::RealizedSamples) = RealizedSamples(copy(r.parent))

reset!(rv::RealizedSamples) = RealizedSamples(reset!(rv.parent))

Base.IteratorEltype(mc::RealizedSamples) = Base.HasEltype()

# function Base.eltype(::Type{Rz}) where {R,S,T,Rz <: RealizedSamples{R,S,T}}
#     eltype(T)
# end

Base.length(r::RealizedSamples) = length(r.parent)

# Base.size(r::RealizedSamples) = ...

Base.IteratorSize(r::RealizedSamples) = Base.IteratorSize(r.parent)

function Base.rand(
    rng::AbstractRNG,
    ::Type{T},
    d::ProductMeasure{G},
) where {T,G<:Base.Generator}
    RealizedSamples(ResettableRNG(rng), marginals(d))
end

function Base.rand(
    rng::AbstractRNG,
    ::Type,
    d::For{T,F,I},
) where {N,T,F,I<:NTuple{N,<:Base.Generator}}
    RealizedSamples(ResettableRNG(rng), marginals(d))
end
