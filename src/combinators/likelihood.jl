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
"""
struct LogLikelihood{F,X}
    f::F
    x::X
end

function show(io::IO, ℓ::LogLikelihood{Tuple{M,NT}}) where {M<:Type, NT <: NamedTuple}
    println(io, "LogLikelihood(",ℓ.f..., ℓ.x, ")")
end

function show(io::IO, ℓ::LogLikelihood)
    println(io, "LogLikelihood(",ℓ.f, ℓ.x, ")")
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
