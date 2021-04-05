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



struct Realized{R,S,T} <: DynamicIterators.DynamicIterator
    seed::R
    rng::S
    iter::T
end
Base.IteratorSize(::Realized) = IsInfinite()
Base.iterate(E::Realized) = dyniterate(E, nothing)
Base.iterate(E::Realized, value) = dyniterate(E, value)

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


# using Soss

# hmm = @model begin
#     ε ~ Exponential() #  transition
#     σ ~ Exponential() # Observation noise
#     x ~ Chain(Normal()) do xj
#         Normal(xj, ε)
#     end

#     y ~ For(x) do xj
#         Normal(xj, σ)
#     end
# end

# using Soss

# mbind = @model μ,κ begin
#     x ~ μ
#     y ~ κ(x)
#     return y
# end

# ⋅(μ,κ) = mbind(μ,κ)

# d =  Cauchy() ⋅ (x -> Normal(μ=x)) ⋅ (x -> Normal(μ=x)) ⋅ (x -> Normal(μ=x))

# rand(d)
# t = xform(d)
# t(randn(4))
# simulate(d)

# # julia> d =  Cauchy() ⋅ (x -> Normal(μ=x)) ⋅ (x -> Normal(μ=x)) ⋅ (x -> Normal(μ=x))
# # ConditionalModel given
# #     arguments    (:μ, :κ)
# #     observations ()
# # @model (μ, κ) begin
# #         x ~ μ
# #         y ~ κ(x)
# #         return y
# #     end



# # julia> rand(d)
# # -3.0414465047589037

# # julia> t = xform(d)
# # TransformVariables.TransformTuple{NamedTuple{(:x, :y), Tuple{TransformVariables.TransformTuple{NamedTuple{(:x, :y), Tuple{TransformVariables.TransformTuple{NamedTuple{(:x, :y), Tuple{TransformVariables.Identity, TransformVariables.Identity}}}, TransformVariables.Identity}}}, TransformVariables.Identity}}}((x = TransformVariables.TransformTuple{NamedTuple{(:x, :y), Tuple{TransformVariables.TransformTuple{NamedTuple{(:x, :y), Tuple{TransformVariables.Identity, TransformVariables.Identity}}}, TransformVariables.Identity}}}((x = TransformVariables.TransformTuple{NamedTuple{(:x, :y), Tuple{TransformVariables.Identity, TransformVariables.Identity}}}((x = asℝ, y = asℝ), 2), y = asℝ), 3), y = asℝ), 4)

# # julia> t(randn(4))
# # (x = (x = (x = -0.24259286698966315, y = 0.278190893626807), y = -1.361907586870645), y = 0.05914265096096323)

# # julia> simulate(d)
# # (value = 5.928939554009484, trace = (x = (value = 5.307006358072237, trace = (x = (value = 3.2023770380851797, trace = (x = 3.677550124255551, y = 3.2023770380851797)), y = 5.307006358072237)), y = 5.928939554009484))
