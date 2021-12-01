using DynamicIterators
import DynamicIterators: dub, dyniterate, evolve
using Base.Iterators: SizeUnknown, IsInfinite

import MeasureBase: For

export Chain

struct Chain{K,M} <: AbstractMeasure
    κ::K
    μ::M
end

Pretty.quoteof(c::Chain) = :(Chain($(Pretty.quoteof(c.κ)), $(Pretty.quoteof(c.μ))))

Base.length(::Chain) = ∞

@inline function basemeasure(mc::Chain)
    Chain(basemeasure ∘ mc.κ, basemeasure(mc.μ))
end

Base.IteratorEltype(mc::Chain) = Base.HasEltype()

Base.eltype(::Type{C}) where {K,M,C<:Chain{K,M}} = eltype(M)

@inline function logdensity_def(mc::Chain, x)
    μ = mc.μ
    ℓ = 0.0
    for xj in x
        ℓ += logdensity_def(μ, xj)
        μ = mc.κ(xj)
    end
    return ℓ
end

DynamicIterators.evolve(mc::Chain, μ) = μ ⋅ mc.κ
DynamicIterators.evolve(mc::Chain) = mc.μ

dyniterate(E::Chain, value) = dub(evolve(E, value))
dyniterate(E::Chain, ::Nothing) = dub(evolve(E))
Base.iterate(E::Chain) = dyniterate(E, nothing)
Base.iterate(E::Chain, value) = dyniterate(E, value)

function DynamicIterators.dyniterate(r::Chain, (x, rng)::Sample)
    μ = r.κ(x)
    y = rand(rng, μ)
    return (μ ↝ y), Sample(y, rng)
end

Base.IteratorSize(::Chain) = IsInfinite()
Base.IteratorSize(::Type{Chain}) = IsInfinite()

function Base.rand(rng::AbstractRNG, T::Type, chain::Chain)
    r = ResettableRNG(rng)
    return RealizedSamples(r, chain)
end

###############################################################################
# DynamicFor

# A `DynamicFor` is produced when `For` is called on a `DynamicIterator`.

struct DynamicFor{T,K,S} <: AbstractMeasure
    κ::K
    iter::S
end

Pretty.quoteof(r::DynamicFor) =
    :(DynamicFor($(Pretty.quoteof(r.κ)), $(Pretty.quoteof(r.iter))))

function DynamicFor(κ::K, iter::S) where {K,S}
    T = typeof(κ(first(iter)))
    DynamicFor{T,K,S}(κ, iter)
end

function Base.rand(rng::AbstractRNG, T::Type, df::DynamicFor)
    r = ResettableRNG(rng)
    return RealizedSamples(r, df)
end

@inline function logdensity_def(df::DynamicFor, y)
    ℓ = 0.0
    for (xj, yj) in zip(df.iter, y)
        ℓ += logdensity_def(df.κ(xj), yj)
    end
    return ℓ
end

Base.eltype(::Type{D}) where {T,D<:DynamicFor{T}} = eltype(T)

Base.IteratorEltype(d::DynamicFor) = Base.HasEltype()

Base.IteratorSize(d::DynamicFor) = Base.IteratorSize(d.iter)

function Base.iterate(d::DynamicFor)
    (x, s) = iterate(d.iter)
    (d.κ(x), s)
end

function Base.iterate(d::DynamicFor, s)
    (x, s) = iterate(d.iter, s)
    (d.κ(x), s)
end

Base.length(d::DynamicFor) = length(d.iter)

For(f, r::Realized) = DynamicFor(f, copy(r))

function Base.rand(rng::AbstractRNG, dfor::DynamicFor)
    r = ResettableRNG(rng)
    return RealizedSamples(r, dfor)
end

function dyniterate(df::DynamicFor, st, args...)
    (val, state) = dyniterate(df.iter, st, args...)
    return (df.κ(val), state)
end

For(f, it::DynamicIterator) = DynamicFor(f, it)

For(f, it::DynamicFor) = DynamicFor(f, it)

function dyniterate(df::DynamicFor, state)
    ϕ = dyniterate(df.iter, state)
    ϕ === nothing && return nothing
    u, state = ϕ
    df.f(u), state
end

function Base.collect(r::Realized)
    next = iterate(r)
    isnothing(next) && return []
    (x, s) = next
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

function testvalue(mc::Chain)
    μ = mc.μ
    κ = mc.κ
    rand(Chain(Dirac ∘ testvalue ∘ κ, (Dirac ∘ testvalue)(μ)))
end
