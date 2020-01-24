module Measures

using MLStyle
using CanonicalTraits
import Distributions
using Reexport

export Measure

abstract type Measure{X} end

@trait IsMeasure{M,X} where {X = eltype(M)}

@implement IsMeasure{M, X} where {X, M <: Measure{X}}


# # This lets us write e.g.
# # @measure Normal
# # making it easier to declare a new measure.
export @measure

"""
    @measure dist

Create a new measure named `dist`. For example, `@measure Normal` is equivalent to

    struct Normal{P, X} <: Measure{X}
        par::P
        Normal(nt::NamedTuple) = new{typeof(nt), _domain(Normal, typeof(nt))}(nt)
    end

    Normal(; kwargs...) = Normal((;kwargs...))
"""
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


# export baseMeasure

# # export ≪

# # @trait DominatedBy{M1 <: Measure{X}, M2 <: Measure{X}} begin
# #     (≪) :: [M1, M2] => Bool
# #     ≪(m1, m2) = false
# # end

export HasDensity

export logdensity

@trait IsMeasure{M, X} >: HasDensity{M, X} where {X = eltype(M)} begin
    logdensity :: [M, X] => Real
end

export measure

@trait IsMeasure{M,X} >: IsMeasurable{M,S,X} where {X = eltype(S)} begin
    measure :: [M, S] => Real
end

export HasRand

@trait IsMeasure{M,X} >: HasRand{M,X} where {X = eltype(M)} begin
    rand :: [M] => eltype(M)
end

# @trait Measure{M,X} >: FiniteMeasure{M,X} where {X = eltype{M}} begin

# end


include("basemeasures/lebesgue.jl")
include("combinators/scale.jl")
# include("combinators/superpose.jl")
include("distributions.jl")

include("probability/normal.jl")

end # module
