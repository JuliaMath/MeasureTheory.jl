export LogLikelihood

@concrete terse struct LogLikelihood{T,X}
    x::X
end

LogLikelihood(T::Type, x::X) where {X} = LogLikelihood{T,X}(x)

LogLikelihood(μ::T, x::X) where {X, T<:AbstractMeasure} = LogLikelihood{T,X}(x)

logdensity(ℓ::LogLikelihood, p) = ℓ(p)

(ℓ::LogLikelihood{T,X})(p) where {T,X} = logdensity(T(p), ℓ.x)

(ℓ::LogLikelihood)(;kwargs...) = ℓ((;kwargs...))
