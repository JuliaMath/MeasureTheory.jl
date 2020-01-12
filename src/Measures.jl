module Measures

using MLStyle
using CanonicalTraits
import Distributions

abstract type Measure{X} end

export domain
domain(::Measure{X}) where {X} = X
domain(::Type{Measure{X}}) where {X} = X



# This lets us write e.g.
# @measure Normal
# making it easier to declare a new measure.
export @measure
macro measure(d)
    d = esc(d)
    return quote
        struct $d{P,X} <: Measure{X}
            par :: P
            $d(nt::NamedTuple) = new{typeof(nt), _domain($d, typeof(nt))}(nt)
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


include("basemeasures/lebesgue.jl")

include("probability/normal.jl")

end # module
