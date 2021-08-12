export Likelihood

@doc raw"""
    Likelihood(M<:ParameterizedMeasure, x)

"Observe" a value `x`, yielding a function from the parameters to ℝ.

Likelihoods are most commonly used in conjunction with an existing _prior_
measure to yield a new measure, the _posterior_. In Bayes's Law, we have

``P(θ|x) ∝ P(θ) P(x|θ)``

Here ``P(θ)`` is the prior. If we consider ``P(x|θ)`` as a function on ``θ``,
then it is called a likelihood.

Since measures are most commonly manipulated using `density` and `logdensity`,
it's awkward to commit a (log-)likelihood to using one or the other. To evaluate
a `Likelihood`, we therefore use `density` or `logdensity`, depending on the
circumstances. In the latter case, it is of course acting as a log-density.

For example,

    julia> ℓ = Likelihood(Normal{(:μ,)}, 2.0)
    Likelihood(Normal{(:μ,), T} where T, 2.0)

    julia> density(ℓ, (μ=2.0,))
    1.0

    julia> logdensity(ℓ, (μ=2.0,))
    -0.0

If, as above, the measure includes the parameter information, we can optionally
leave it out of the second argument in the call to `density` or `logdensity`. 

    julia> density(ℓ, 2.0)
    1.0

    julia> logdensity(ℓ, 2.0)
    -0.0

With several parameters, things work as expected:
    
    julia> ℓ = Likelihood(Normal{(:μ,:σ)}, 2.0)
    Likelihood(Normal{(:μ, :σ), T} where T, 2.0)
    
    julia> logdensity(ℓ, (μ=2, σ=3))
    -1.0986122886681098
    
    julia> logdensity(ℓ, (2,3))
    -1.0986122886681098
    
    julia> logdensity(ℓ, [2, 3])
    -1.0986122886681098

---------

    Likelihood(M<:ParameterizedMeasure, constraint::NamedTuple, x)

In some cases the measure might have several parameters, and we may want the
(log-)likelihood with respect to some subset of them. In this case, we can use
the three-argument form, where the second argument is a constraint. For example,

    julia> ℓ = Likelihood(Normal{(:μ,:σ)}, (σ=3.0,), 2.0)
    Likelihood(Normal{(:μ, :σ), T} where T, (σ = 3.0,), 2.0)

Similarly to the above, we have

    julia> density(ℓ, (μ=2.0,))
    0.3333333333333333

    julia> logdensity(ℓ, (μ=2.0,))
    -1.0986122886681098

    julia> density(ℓ, 2.0)
    0.3333333333333333

    julia> logdensity(ℓ, 2.0)
    -1.0986122886681098

-----------------------

Finally, let's return to the expression for Bayes's Law, 

``P(θ|x) ∝ P(θ) P(x|θ)``

The product on the right side is computed pointwise. To work with this in
MeasureTheory, we have a "pointwise product" `⊙`, which takes a measure and a
likelihood, and returns a new measure, that is, the unnormalized posterior that has density ``P(θ) P(x|θ)`` with respect to the base measure of the prior.

For example, say we have

    μ ~ Normal()
    x ~ Normal(μ,σ)
    σ = 1

and we observe `x=3`. We can compute the posterior measure on `μ` as

    julia> post = Normal() ⊙ Likelihood(Normal{(:μ, :σ)}, (σ=1,), 3)
    Normal() ⊙ Likelihood(Normal{(:μ, :σ), T} where T, (σ = 1,), 3)

    julia> logdensity(post, 2)
    -2.5
"""
struct Likelihood{F,S,X}
    k::Kernel{F,S}
    x::X
end

Likelihood(μ::AbstractMeasure, x) = Likelihood(kernel(μ),x)

Likelihood(::Type{M}, x) where {M<:AbstractMeasure} = Likelihood(kernel(M),x)

function Base.show(io::IO, ℓ::Likelihood) 
    io = IOContext(io, :compact => true)
    k, x = ℓ.k, ℓ.x
    print(io, "Likelihood(",k,", ", x, ")")
end


function logdensity(ℓ::Likelihood, p) 
    return logdensity(ℓ.k(p), ℓ.x)
end
