using Random

struct ResettableRNG{R,S} <: Random.AbstractRNG
    rng::R
    seed::S
end

function Base.show(io::IO, r::ResettableRNG)
    io = IOContext(io, :compact => true)
    print(io, "ResettableRNG(::", constructor(r.rng), ", ", r.seed, ")")
end

function reset!(r::ResettableRNG)
    Random.seed!(r.rng, r.seed)
end


# for f in [
#     :(Base.rand)
#     :(Base.randn)
#     :(Random.rand!)
#     :(Random.randcycle)
#     :(Random.randcycle!)
#     :(Random.randexp)
#     :(Random.randexp!)
#     :(Random.randn!)
#     :(Random.randperm)
#     :(Random.randperm!)
#     :(Random.randstring)
#     :(Random.randsubseq)
#     :(Random.randsubseq!)
#     :(Random.shuffle)
#     :(Random.shuffle!)
#     :(Random.seed!)
# ]    
#     @eval $f(r::ResettableRNG, args...) = $f(r.rng, args...)
# end


# Base.rand(r::ResettableRNG, d::AbstractMeasure) = rand(r.rng, d)
# Base.rand(r::ResettableRNG, ::Type{T}, d::AbstractMeasure) where {T} = rand(r.rng, T, d)
# Base.rand(r::ResettableRNG) = rand(r.rng, Float64)
Base.rand(r::ResettableRNG, ::Type{T}) where {T} =  rand(r.rng, T)
Base.randn(r::ResettableRNG, ::Type{T}) where {T} =  randn(r.rng, T)
Base.randn(r::ResettableRNG, ::Type{Float64}) where {T} =  randn(r.rng, Float64)

function Base.iterate(r::ResettableRNG)
    reset!(r)
    return (rand(r), nothing)
end

Base.iterate(r::ResettableRNG, _) = (rand(r), nothing)
Base.IteratorSize(r::ResettableRNG) = Base.IsInfinite()
