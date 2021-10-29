using ConcreteStructs
using DynamicIterators
using DynamicIterators: dub
using Base.Iterators: SizeUnknown, IsInfinite
import DynamicIterators: dyniterate, evolve

import MeasureBase: For

export Chain

@concrete terse struct Chain{K,M} <: AbstractMeasure
    κ::K
    μ::M
end

Pretty.quoteof(c::Chain) = :(Chain($(Pretty.quoteof(c.κ)), $(Pretty.quoteof(c.μ))))



function basemeasure(mc::Chain)
    Chain(basemeasure ∘ mc.κ, basemeasure(mc.μ))
end

Base.IteratorEltype(mc::Chain) = Base.HasEltype()

Base.eltype(::Type{C}) where {K,M,C<:Chain{K,M}} = eltype(M)

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
    rng = deepcopy(rng)
    u = rand(rng, μ)
    return u, Sample(u, rng)
end
Base.IteratorSize(::Chain) = IsInfinite()
Base.IteratorSize(::Type{Chain}) = IsInfinite()


@concrete terse struct Realized{R,S,T} <: DynamicIterators.DynamicIterator
    rng::ResettableRNG{R,S}
    iter::T
end

Base.show(io::IO, r::Realized) = Pretty.pprint(io, r)

Pretty.quoteof(r::Realized) = :(Realized($(Pretty.quoteof(r.rng)), $(Pretty.quoteof(r.iter))))

Base.IteratorEltype(mc::Realized) = Base.HasEltype()

function Base.eltype(::Type{Rz}) where {R,S,T,Rz <: Realized{R,S,T}}
    eltype(T)
end

Base.length(r::Realized) = length(r.iter)

Base.size(r::Realized) = size(r.iter)

Base.IteratorSize(::Type{Rz}) where {R,S,T, Rz <: Realized{R,S,T}} = Base.IteratorSize(T)
Base.IteratorSize(r::Rz) where {R,S,T, Rz <: Realized{R,S,T}} = Base.IteratorSize(r.iter)


function Base.iterate(rv::Realized{R,S,T}) where {R,S,T}
    if static_hasmethod(evolve, Tuple{T})
        dyniterate(rv, nothing)
    else
        !isnothing(rv.rng.seed) && reset!(rv.rng)
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
    !isnothing(rv.rng.seed) && reset!(rv.rng)
    μ = evolve(rv.iter)
    x = rand(rv.rng, μ)
    x, Sample(x, rv.rng)
end
function dyniterate(rv::Realized, u::Sample)
    dyniterate(rv.iter, u)
end

function Base.rand(rng::AbstractRNG, T::Type, chain::Chain)
    seed = rand(rng, UInt)
    r = ResettableRNG(rng, seed)
    return Realized(r, chain)
end

###############################################################################
# DynamicFor

# A `DynamicFor` is produced when `For` is called on a `DynamicIterator`.

@concrete terse struct DynamicFor{T,K,S} <: AbstractMeasure
    κ ::K
    sampler :: S        
end

Pretty.quoteof(r::DynamicFor) = :(DynamicFor($(Pretty.quoteof(r.κ)), $(Pretty.quoteof(r.sampler))))

function DynamicFor(κ::K,sampler::S) where {K,S}
    T = typeof(κ(first(sampler)))
    DynamicFor{T,K,S}(κ,sampler)
end

function Base.rand(rng::AbstractRNG, T::Type, df::DynamicFor)
    rng = deepcopy(rng)
    seed = rand(rng, UInt)
    r = ResettableRNG(rng, seed)
    return Realized(r, df)
end

function logdensity(df::DynamicFor, y)
    ℓ = 0.0
    for (xj, yj) in zip(df.sampler, y)
        ℓ += logdensity(df.κ(xj), yj)
    end
    return ℓ
end

Base.eltype(::Type{D}) where {T,D<:DynamicFor{T}} = eltype(T)

Base.IteratorEltype(d::DynamicFor) = Base.HasEltype()

Base.IteratorSize(d::DynamicFor) = Base.IteratorSize(d.sampler)

function Base.iterate(d::DynamicFor)
    (x,s) = iterate(d.sampler)
    (d.κ(x), s)
end

function Base.iterate(d::DynamicFor, s)
    (x,s) = iterate(d.sampler, s)
    (d.κ(x), s)
end

Base.length(d::DynamicFor) = length(d.sampler)


For(f, r::Realized) = DynamicFor(f,r)

function Base.rand(rng::AbstractRNG, dfor::DynamicFor)
    seed = rand(rng, UInt)
    rng = deepcopy(rng)
    r = ResettableRNG(rng, seed)
    return Realized(r, dfor)
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

function testvalue(mc::Chain)
    μ = mc.μ
    κ = mc.κ
    rand(Chain(Dirac ∘ testvalue ∘ κ, (Dirac ∘ testvalue)(μ)))
end
