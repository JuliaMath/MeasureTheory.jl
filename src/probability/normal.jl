
# Normal distribution


export Normal

@measure Normal
baseMeasure(::Type{Normal{P,X}}) where {P,X} = Lebesgue(X)

Normal(μ::Real, σ::Real) = Normal(μ=μ, σ=σ)


# Standard normal

_domain(::Type{Normal}, ::Type{NamedTuple{(),Tuple{}}}) = Float64


@implement Density{Normal{P,X},X} where {X <: Real, P <: NamedTuple{(),Tuple{}}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π) + x^2) / 2  
end

# μ, σ


_domain(::Type{Normal}, ::Type{NamedTuple{(:μ, :σ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)


@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π)) / 2 - log(d.par.σ)  - (x - d.par.μ)^2 / (2 * d.par.σ^2)
end

# μ, τ

_domain(::Type{Normal}, ::Type{NamedTuple{(:μ, :τ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)

@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :τ)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π) - log(d.par.τ)  + d.par.τ * (x - d.par.μ)^2) / 2
end


# σ


_domain(::Type{Normal}, ::Type{NamedTuple{(:σ,), Tuple{B}}}) where {B} = B


@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:σ,)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π)) / 2 - log(d.par.σ)  - x^2 / (2 * d.par.σ^2)
end

# τ

_domain(::Type{Normal}, ::Type{NamedTuple{(:τ,), Tuple{B}}}) where {B} = B

@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:τ,)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π) - log(d.par.τ)  + d.par.τ * x^2) / 2
end

# logdensity(Normal(σ=0.5), 4.0) 
#  == logdensity(Normal(μ=0.0, σ=0.5), 4.0) 
#  == logdensity(Normal(μ=0.0, τ=4.0), 4.0) 
#  == logdensity(Normal(τ=4.0), 4.0)

