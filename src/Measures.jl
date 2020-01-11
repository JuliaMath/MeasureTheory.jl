module Measures

using MLStyle
using CanonicalTraits
import Distributions

abstract type Measure{X} end

export domain
domain(::Measure{X}) where {P,X} = X
domain(::Type{Measure{X}}) where {P,X} = X



# This lets us write e.g.
# @measure Normal
# making it easier to declare a new measure.
# We should extend this to allow
# @measure Normal(μ,σ)
# which could add a method
# Normal(μ,σ) = Normal(μ=μ,σ=σ)
# for default names
export @measure
macro measure(d)
    d = esc(d)
    return quote
        struct $d{P,X} <: Measure{X}
            par :: P
            $d(nt::NamedTuple) = new{typeof(nt), getX(Normal, typeof(nt))}(nt)
        end
        $d(;kwargs...) = $d((;kwargs...))
    end    
end


export baseMeasure
export logdensity

@trait Density{M,X} where {X = domain{M}} begin
    baseMeasure :: [M] => Measure{X}
    logdensity :: [M, X] => Real
end

@trait Measurable{M,S} begin
    measure :: [M, S] => Real
end

# Lebesgue measure

export Lebesgue
struct Lebesgue{X} <: Measure{X} end
Lebesgue(X) = Lebesgue{X}()

# Normal distribution

export Normal

@measure Normal


getX(::Type{Normal}, ::Type{NamedTuple{(:μ, :σ),Tuple{A,B}}}) where {A,B} = promote_type(A,B)

Normal(μ::Real, σ::Real) = Normal(μ=μ, σ=σ)

baseMeasure(::Type{Normal{P,X}}) where {P,X} = Lebesgue(X)



@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = -log(d.par.σ) - log(2)/2 - log(π)/2 - (x - d.par.μ)^2 / (2 * d.par.σ^2)
end

end # module
