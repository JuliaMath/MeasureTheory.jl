using ConcreteStructs
using DynamicIterators
using DynamicIterators: dub
using Base.Iterators: SizeUnknown, IsInfinite
import DynamicIterators: dyniterate, evolve

export Chain

@concrete terse struct Chain <: AbstractMeasure
    κ
    μ
end

function logdensity(mc::Chain, x)
    μ = mc.μ
    ℓ = 0.0
    for xj in x
        ℓ += logdensity(μ, xj)
        μ = mc.κ(xj)
    end
    return ℓ
end

DynamicIterators.evolve(mc::Chain, μ) =  μ ⋅ mc.κ
DynamicIterators.evolve(mc::Chain) =  mc.μ

dyniterate(E::Chain, value) = dub(evolve(E, value))
dyniterate(E::Chain, ::Nothing) = dub(evolve(E))
Base.iterate(E::Chain) = dyniterate(E, nothing)
Base.iterate(E::Chain, value) = dyniterate(E, value)

function DynamicIterators.dyniterate(r::Chain, (u,rng)::Sample)
    μ = r.κ(u) 
    u = rand(rng, μ)
    return u, Sample(u, rng)
end
Base.IteratorSize(::Chain) = IsInfinite()
Base.IteratorSize(::Type{Chain}) = IsInfinite()


struct Realized{R,S,T} <: DynamicIterators.DynamicIterator
    seed::R
    rng::S
    iter::T
end

Base.size(r::Realized) = size(r.iter)

Base.IteratorSize(::Type{Rz}) where {R,S,T, Rz <: Realized{R,S,T}} = Base.IteratorSize(T)
Base.IteratorSize(r::Rz) where {R,S,T, Rz <: Realized{R,S,T}} = Base.IteratorSize(r.iter)


function Base.iterate(rv::Realized{R,S,T}) where {R,S,T}
    if static_hasmethod(evolve, Tuple{T})
        dyniterate(rv, nothing)
    else
        !isnothing(rv.seed) && Random.seed!(rv.rng, rv.seed)
        μ,s = iterate(rv.iter)
        x = rand(rv.rng, μ)
        x,s
    end
end


function Base.iterate(rv::Realized{R,S,T}, s) where {R,S,T}
    if static_hasmethod(evolve, Tuple{T})
        dyniterate(rv, s)
    else
        μs = iterate(rv.iter, s)
        isnothing(μs) && return nothing
        (μ,s) = μs
        x = rand(rv.rng, μ)
        return x,s
    end
end


function dyniterate(rv::Realized, ::Nothing)
    !isnothing(rv.seed) && Random.seed!(rv.rng, rv.seed)
    μ = evolve(rv.iter)
    x = rand(rv.rng, μ)
    x, Sample(x, rv.rng)
end
function dyniterate(rv::Realized, u::Sample)
    dyniterate(rv.iter, u)
end

function Base.rand(rng::AbstractRNG, T::Type, chain::Chain)
    seed = rand(rng, UInt)
    return Realized(seed, copy(rng), chain)
end

###############################################################################
# DynamicFor

# A `DynamicFor` is produced when `For` is called on a `DynamicIterator`.

@concrete terse struct DynamicFor <: AbstractMeasure
    κ
    sampler
end

For(f, r::Realized) = DynamicFor(f,r)

function Base.rand(rng::AbstractRNG, dfor::DynamicFor)
    seed = rand(rng, UInt)
    return Realized(seed, copy(rng), dfor)
end

function dyniterate(df::DynamicFor, st, args...)
    (val, state) = dyniterate(df.iter, st, args...)
    return (df.κ(val), state)
end

For(f, it::DynamicIterator) = DynamicFor(f, it)

For(f, it::DynamicFor) = DynamicFor(f, it)

function dyniterate(fr::DynamicFor, state)
      ϕ = dyniterate(fr.iter, state)
      ϕ === nothing && return nothing
      u, state = ϕ
      fr.f(u), state
end

function Base.collect(r::Realized)
    next = iterate(r)
    isnothing(next) && return []
    (x,s) = next
    a = similar(r.iter, typeof(x))

    i = 1
    @inbounds a[i] = x
    while !isnothing(next)
        (x, s) = next
        @inbounds a[i] = x
        i += 1
        next = iterate(r, s)
    end
    return a
end 
