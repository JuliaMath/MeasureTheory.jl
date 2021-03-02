export Likelihood

struct Likelihood{T,X} <: AbstractMeasure
    x::X
end

Likelihood(T::Type, x::X) where {X} = Likelihood{T,X}(x)

function logdensity(μ::Likelihood{T,X}, p::NamedTuple) where {T,X}
    logdensity(T(p), μ.x)
end
