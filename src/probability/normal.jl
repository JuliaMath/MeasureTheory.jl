
# Normal distribution

import StatsFuns
using Intervals
export Normal

@measure Normal
baseMeasure(::Normal{P,X}) where {P,X} = Lebesgue(X)

Normal(μ::Real, σ::Real) = Normal(μ=μ, σ=σ)


# # Standard normal

_domain(::Type{Normal}, ::Type{NamedTuple{(),Tuple{}}}) = Float64


# # μ, σ


_domain(::Type{Normal}, ::Type{NamedTuple{(:μ, :σ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)

function logdensity(d::Normal{P,X} , x::X) where {P <: NamedTuple{(:μ, :σ)}, X <: Real}
    return - (log(2) + log(π)) / 2 - log(d.par.σ)  - (x - d.par.μ)^2 / (2 * d.par.σ^2)
end


 
# @implement HasDensity{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
#     baseMeasure(d) = Lebesgue(X)
#     logdensity(d, x) = - (log(2) + log(π)) / 2 - log(d.par.σ)  - (x - d.par.μ)^2 / (2 * d.par.σ^2)
# end

# # μ, τ

# _domain(::Type{Normal}, ::Type{NamedTuple{(:μ, :τ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)

# @implement HasDensity{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :τ)}} begin
#     baseMeasure(d) = Lebesgue(X)
#     logdensity(d, x) = - (log(2) + log(π) - log(d.par.τ)  + d.par.τ * (x - d.par.μ)^2) / 2
# end


# # σ


# _domain(::Type{Normal}, ::Type{NamedTuple{(:σ,), Tuple{B}}}) where {B} = B


# @implement HasDensity{Normal{P,X},X} where {X, P <: NamedTuple{(:σ,)}} begin
#     baseMeasure(d) = Lebesgue(X)
#     logdensity(d, x) = - (log(2) + log(π)) / 2 - log(d.par.σ)  - x^2 / (2 * d.par.σ^2)
# end

# # τ

# _domain(::Type{Normal}, ::Type{NamedTuple{(:τ,), Tuple{B}}}) where {B} = B

# @implement HasDensity{Normal{P,X},X} where {X, P <: NamedTuple{(:τ,)}} begin
#     baseMeasure(d) = Lebesgue(X)
#     logdensity(d, x) = - (log(2) + log(π) - log(d.par.τ)  + d.par.τ * x^2) / 2
# end

# @implement HasRand{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
#     rand(d) = Distributions.rand(Normal(d.par.μ,d.par.σ))
# end


# @implement IsMeasurable{Normal{P,X},S,X} where {X <: Real, S <: Interval{X}, P <: NamedTuple{(:μ, :σ)}} begin
#     measure(d, i) = StatsFuns.normcdf(d.par.μ, d.par.σ, i.last) - StatsFuns.normcdf(d.par.μ, d.par.σ, i.first)
# end

# logdensity(Normal(σ=0.5), 4.0) 
#  == logdensity(Normal(μ=0.0, σ=0.5), 4.0) 
#  == logdensity(Normal(μ=0.0, τ=4.0), 4.0) 
#  == logdensity(Normal(τ=4.0), 4.0)
