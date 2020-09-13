module Measures

using MLStyle
import Distributions
using Reexport
using BangBang
using SimplePosets

# include("traits.jl")

export Measure

abstract type Measure{X} end

@traitdef IsMeasure{X}

@traitimpl IsMeasure{Distribution}

@traitfn

Base.eltype(μ::Measure{X}) where {X} = X

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

include("basemeasures.jl")
include("basemeasures/lebesgue.jl")
include("combinators/scale.jl")
include("distributions.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
# include("probability/normal.jl")



end # module
