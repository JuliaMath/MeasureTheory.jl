using ConcreteStructs
using DynamicIterators

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
DynamicIterators.evolve(mc::Chain, ::Nothing) =  mc.μ

@concrete terse struct SampleInit
    rng
    seed
end

function dyniterate(mc::Chain, (u,)::Start{<:SampleInit})
    !isnothing(u.seed) && Random.seed!(u.rng, u.seed)
    x = rand(u.rng, mc.μ)
    1=>x, Sample(1=>x, u.rng)
end

function dyniterate(mc::Chain, (i=>x,rng)::Sample)
    xnew = rand(rng, mc.κ(x))
    xnew, Sample(i+1 => xnew, rng)
end

Base.IteratorSize(::Chain) = Base.IsInfinite()

function Base.rand(rng::AbstractRNG, mc::Chain)
    rng = copy(rng)
    seed = rand(rng, UInt)
    return dyniterate(mc, Start(SampleInit(rng, seed)))
end

const DynamicSampler{T,S,R} = Tuple{T, DynamicIterators.Sample{S,R}}


###############################################################################
# DynamicFor

# A `DynamicFor` is produced when `For` is called on a `DynamicIterator`.

@concrete terse struct DynamicFor{F,T,S,R} <: AbstractMeasure
    f :: F
    sampler::DynamicSampler{T,S,R}
end

For(f, s::DynamicSampler) = DynamicFor(f,s)

function dyniterate(df::DynamicFor, st, args...)
    (val, state) = dyniterate(df.iter, st, args...)
    return (f(val), state)
end

For(f, it::DynamicIterator) = DynamicFor(f, it)

For(f, it::DynamicFor) = DynamicFor(f, it)

function dyniterate(fr::DynamicFor, state)
      ϕ = dyniterate(fr.iter, state)
      ϕ === nothing && return nothing
      u, state = ϕ
      fr.f(u), state
end


# @concrete terse struct RandChain <: Evolution
#     rng
#     mc
# end

# Base.IteratorSize(::RandChain) = Base.IsInfinite()

# function DynamicIterators.evolve(rc::RandChain, x) 
#     rng = rc.rng
#     rand(rng, rc.mc.κ(x))
# end

# function DynamicIterators.evolve(rc::RandChain, ::Nothing)
#     rng = rc.rng
#     μ0 = rc.mc.μ
#     rand(rng, μ0)
# end


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
