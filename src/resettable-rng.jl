using Random
export ResettableRNG

export reset!
struct ResettableRNG{R,S} <: Random.AbstractRNG
    rng::R
    seed::S

    function ResettableRNG(rng, seed::S) where {S}
        rng = copy(rng)
        R2 = typeof(rng) 
        new{R2,S}(rng, seed)
    end

    function ResettableRNG(rng)
        seed = rand(rng, UInt)
        S = typeof(seed)
        rng = copy(rng)
        R2 = typeof(rng) 
        new{R2,UInt}(copy(rng), seed)
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
    # @info "Calling reset! on $r"
    Random.seed!(r.rng, r.seed)
end


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

Base.iterate(r::ResettableRNG) = iterate(r, nothing)

function Base.iterate(r::ResettableRNG, ::Nothing)
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
