export LogLikelihood

"""
    LogLikelihood(M<:ParameterizedMeasure, x)

"Observe" a value `x`, yielding a function from the parameters to ℝ.

EXAMPLE:
    julia> LogLikelihood(Normal, 3)(μ=2)
    -0.5

If `M` has a name, we don't need one in the argument.

    julia> LogLikelihood(Normal{(:μ,)}, 3)(2)
    -0.5

---------

    LogLikelihood(M<:ParameterizedMeasure, constraint::NamedTuple, x)

As above, but merge the new values with `constraint` before evaluation

    julia> LogLikelihood(Normal, (μ=2,), 5)((σ=1,))
    -4.5

We can again provide named in the first argument:

    julia> LogLikelihood(Normal{(:μ,:σ)}, (μ=2,), 5)(1)
    -4.5

--------

This last form is especially useful for probabilistic programming. Suppose we have a
`μ ~ Normal()` and `x ~ Normal(μ,σ)`. If we observe `x=3` and, say, `σ=1`, then
the posterior density at `μ` is

    post = Normal() ⊙ LogLikelihood(Normal{(:μ,:σ)}, (σ=1,), 3)

or, for example,

    julia> logdensity(post, 2)
    -2.5
"""
struct LogLikelihood{F,X}
    f::F
    x::X
end

function Base.show(io::IO, ℓ::LogLikelihood{Tuple{M,NamedTuple{N,T}}}) where {M<:Type, N,T}
    (m,c) = ℓ.f
    println(io, "LogLikelihood(",m,", ",c,", ", ℓ.x, ")")
end

function Base.show(io::IO, ℓ::LogLikelihood)
    println(io, "LogLikelihood(",ℓ.f, ", ", ℓ.x, ")")
end

LogLikelihood(μ::M, x) where {M<:AbstractMeasure} = LogLikelihood(M, x)

function LogLikelihood(::Type{M}, constraint::NamedTuple, x) where {M <: ParameterizedMeasure}
    LogLikelihood((M, constraint), x)
end

function LogLikelihood(μ::M, constraint::NamedTuple, x) where {M<:AbstractMeasure}
    LogLikelihood((M, constraint), x)
end

function (ℓ::LogLikelihood{Tuple{M,NT}})(p::NamedTuple) where {M<:Type, NT <: NamedTuple}
    (D, constraint) = ℓ.f
    return logdensity(D(merge(p, constraint)), ℓ.x)
end

function (ℓ::LogLikelihood{Tuple{M,NT}})(p) where {M<:Type, NT <: NamedTuple}
    freevar = params(ℓ.f...)
    ℓ(NamedTuple{freevar}(p))
end

logdensity(ℓ::LogLikelihood, p) = ℓ(p)

(ℓ::LogLikelihood)(p) = logdensity(ℓ.f(p), ℓ.x)

(ℓ::LogLikelihood)(;kwargs...) = ℓ((;kwargs...))
