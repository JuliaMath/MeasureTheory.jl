export Likelihood

@concrete terse struct Likelihood{T,X}
    x::X
end

Likelihood(T::Type, x::X) where {X} = Likelihood{T,X}(x)

Likelihood(μ::T, x::X) where {X, T<:AbstractMeasure} = Likelihood{T,X}(x)

logdensity(ℓ::Likelihood, p) = ℓ(p)

(ℓ::Likelihood{T,X})(p) where {T,X} = logdensity(T(p), ℓ.x)

(ℓ::Likelihood)(;kwargs...) = ℓ((;kwargs...))
