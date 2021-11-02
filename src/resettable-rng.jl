using Random
export ResettableRNG

export reset!
struct ResettableRNG{R,S} <: Random.AbstractRNG
    rng::R
    seed::S

    function ResettableRNG(rng::R, seed::S) where {R,S}
        new{R,S}(copy(rng), seed)
    end

    function ResettableRNG(rng::R) where {R}
        seed = rand(rng, UInt)
        S = typeof(seed)
        new{R,UInt}(copy(rng), seed)
    end
end

function Base.copy(r::ResettableRNG)
    ResettableRNG(r.rng, copy(r.seed))
end

function Base.show(io::IO, r::ResettableRNG)
    io = IOContext(io, :compact => true)
    print(io, "ResettableRNG(::", constructor(r.rng), ", ", r.seed, ")")
end

function reset!(r::ResettableRNG)
    @info "Calling reset! on $r"
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

import Random
using Random: randexp
using InteractiveUtils: subtypes

for T in vcat(subtypes(Signed), subtypes(Unsigned), subtypes(AbstractFloat))
    isbitstype(T) || continue

    @eval begin
        
        function Base.rand(r::ResettableRNG, ::Type{$T}) 
            rand(r.rng, $T)
        end

        function Base.randn(r::ResettableRNG, ::Type{$T}) 
            randn(r.rng, $T)
        end

        function Random.randexp(r::ResettableRNG, ::Type{$T})
            randexp(r.rng, $T)
        end
    end
end

function Base.iterate(r::ResettableRNG)
    @info "Calling `iterate(::ResettableRNG)"
    r = copy(r)
    reset!(r)
    return (rand(r), nothing)
end

Base.iterate(r::ResettableRNG, _) = (rand(r), nothing)
Base.IteratorSize(r::ResettableRNG) = Base.IsInfinite()

function Random.Sampler(r::Type{R}, s::Random.Sampler, rep::Random.Repetition) where {R<:ResettableRNG}
    return Random.Sampler(r.rng, s, r)
end

function Base.rand(r::ResettableRNG, sp::Random.Sampler)
    rand(r.rng, sp)
end
