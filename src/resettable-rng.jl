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
    print(io, "ResettableRNG(::", constructorof(typeof(r.rng)), ", ", r.seed, ")")
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

function Random.Sampler(
    r::Type{R},
    s::Random.Sampler,
    rep::Random.Repetition,
) where {R<:ResettableRNG}
    return Random.Sampler(r.rng, s, r)
end

UIntBitsTypes = [UInt128, UInt16, UInt32, UInt64, UInt8]
IntBitsTypes = [Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8]
FloatBitsTypes = [Float16, Float32, Float64]

for I in IntBitsTypes
    for T in [
        Random.SamplerTrivial{Random.UInt104Raw{I}} 
        Random.SamplerTrivial{Random.UInt10Raw{I}} 
    ] 
        @eval begin
            function Base.rand(r::ResettableRNG, sp::$T)
                rand(r.rng, sp)
            end
        end
    end
end

for U in UIntBitsTypes
    for I in IntBitsTypes
        for T in [
            Random.SamplerRangeInt{T,U} where T<:Union{IntBitsTypes...}
            Random.SamplerRangeFast{U,I}
        ]
            @eval begin
                function Base.rand(r::ResettableRNG, sp::$T)
                    rand(r.rng, sp)
                end
            end
        end
    end
end



for T in [
    Random.Sampler
    Random.SamplerBigInt
    Random.SamplerTag{<:Set, <:Random.Sampler}
    # Random.SamplerTrivial{Random.CloseOpen01{T}} where {T<:FloatBitsTypes}
    # Random.SamplerTrivial{Random.UInt23Raw{UInt32}}
    Random.UniformT
    Random.SamplerSimple{T, S, E} where {E, S, T<:Tuple} 
    Random.SamplerType{T} where T<:AbstractChar
    Random.SamplerTrivial{Tuple{A}} where A
    Random.SamplerSimple{Tuple{A, B, C}, S, E} where {E, S, A, B, C}
    Random.SamplerSimple{<:AbstractArray, <:Random.Sampler}
    Random.Masked
    Random.SamplerSimple{BitSet, <:Random.Sampler}
    Random.SamplerTrivial{<:Random.UniformBits{T}, E} where {E,T}
] 
    @eval begin
        function Base.rand(r::ResettableRNG, sp::$T)
            rand(r.rng, sp)
        end
    end
end
