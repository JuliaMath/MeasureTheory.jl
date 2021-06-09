export LogLikelihood

@concrete terse struct LogLikelihood
    f
    x
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
