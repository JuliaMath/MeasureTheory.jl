
# Normal distribution

export Normal

@measure Normal


getX(::Type{Normal}, ::Type{NamedTuple{(:μ, :σ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)

Normal(μ::Real, σ::Real) = Normal(μ=μ, σ=σ)

baseMeasure(::Type{Normal{P,X}}) where {P,X} = Lebesgue(X)



@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π)) / 2 - log(d.par.σ)  - (x - d.par.μ)^2 / (2 * d.par.σ^2)
end