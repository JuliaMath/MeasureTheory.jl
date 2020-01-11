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

# @measure Normal

struct Normal{P,X} <: Measure{X}
    par :: P
end

Normal(;kwargs...) = Normal((;kwargs...))



function Normal(μ::Pμ, σ::Pσ) where {Pμ <: Real, Pσ <: Real}
    let T = promote_type(Float64, Pμ, Pσ)
        NT = NamedTuple{(:μ, :σ),Tuple{T, T}}
        Normal{NT,T}((μ = T(μ), σ = T(σ)))
    end
end

baseMeasure(::Type{Normal{P,X}}) where {P,X} = Lebesgue(X)



@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = -log(d.par.σ) - log(2)/2 - log(π)/2 - (x - d.par.μ)^2 / (2 * d.par.σ^2)
end

end # module
