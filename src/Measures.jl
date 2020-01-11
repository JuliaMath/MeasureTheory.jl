module Measures

using MLStyle
using CanonicalTraits
import Distributions

abstract type Measure{X} end

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
macro measure(expr)
    @match expr begin
        d => begin
            d = esc(d)
            return quote
                struct $d{P,X} <: Measure{X}
                    par :: P
                end

                $d(;kwargs...) = $d((;kwargs...))
            end
        end
    end
    
end




@trait Density{M,X} where {X = domain{M}} begin
    baseMeasure :: [M] => Measure{X}
    logdensity :: [M, X] => Real
end

@trait Measurable{M,S} begin
    measure :: [M, S] => Real
end

# Lebesgue measure

struct Lebesgue{X} <: Measure{X} end
Lebesgue(X) = Lebesgue{X}()

# Normal distribution

export Normal
@measure Normal

baseMeasure(::Type{Normal{P,X}}) where {P,X} = Lebesgue(X)



@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = -log(d.σ) - log(2)/2 - log(π)/2 - (x - d.μ)^2 / (2 * d.σ^2)
end

end # module
