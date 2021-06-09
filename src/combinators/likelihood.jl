export LogLikelihood

@concrete terse struct LogLikelihood
    f
    x
end

LogLikelihood(μ::M, x) where {M<:AbstractMeasure} = LogLikelihood(constructor(M), x)

function LogLikelihood(M::Type, constraint::NamedTuple, x) 
    LogLikelihood((constructor(M), constraint), x)
end

function LogLikelihood(μ::M, constraint::NamedTuple, x) where {M<:AbstractMeasure}
    LogLikelihood((constructor(M), constraint), x)
end

function (ℓ::LogLikelihood{Tuple{M,NT}})(p) where {M<:Type, NT <: NamedTuple}
    (D, constraint) = ℓ.f
    return logdensity(D(merge(p, constraint)), ℓ.x)
end

logdensity(ℓ::LogLikelihood, p) = ℓ(p)

(ℓ::LogLikelihood)(p) = logdensity(ℓ.f(p), ℓ.x)

(ℓ::LogLikelihood)(;kwargs...) = ℓ((;kwargs...))
