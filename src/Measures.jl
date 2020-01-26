module Measures

using MLStyle
using CanonicalTraits
import Distributions
using Reexport


include("traits.jl")

export Measure

abstract type Measure{X} end



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

include("basemeasures/lebesgue.jl")
include("combinators/scale.jl")
# include("combinators/superpose.jl")
include("distributions.jl")

include("probability/normal.jl")

end # module
