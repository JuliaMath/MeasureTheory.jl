using ConcreteStructs
using DynamicIterators
using DynamicIterators: dub
import Base

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

Base.IteratorSize(::Chain) = Base.IsInfinite()

###############################################################################
# RandChain

Base.rand(rng::AbstractRNG, mc::Chain) = RandChain(rng, mc)

@concrete terse struct RandChain <: Evolution
    rng
    mc
end

Base.IteratorSize(::RandChain) = Base.IsInfinite()

function DynamicIterators.evolve(rc::RandChain, x) 
    rng = rc.rng
    rand(rng, rc.mc.κ(x))
end

function DynamicIterators.evolve(rc::RandChain, ::Nothing)
    rng = rc.rng
    μ0 = rc.mc.μ
    rand(rng, μ0)
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
