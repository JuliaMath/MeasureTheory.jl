module Measures

using MLStyle
import Distributions
using Reexport
using SimplePosets
using SimpleTraits

import Distributions
const Dists = Distributions
# include("traits.jl")



abstract type AbstractMeasure{X} end




@traitdef IsMeasure{X}

# @traitimpl IsMeasure{Dists.Distribution}

Base.eltype(μ::AbstractMeasure{X}) where {X} = X



logdensity(μ::Dists.Distribution, x) = Dists.logpdf(μ,x)

density(μ::Dists.Distribution, x) = Dists.pdf(μ,x)



# # This lets us write e.g.
# # @measure Normal
# # making it easier to declare a new measure.
export @measure

"""
    @measure dist

Create a new measure named `dist`. For example, `@measure Normal` is equivalent to

    struct Normal{P, X} <: AbstractMeasure{X}
        par::P
        Normal(nt::NamedTuple) = new{typeof(nt), _domain(Normal, typeof(nt))}(nt)
    end

    Normal(; kwargs...) = Normal((;kwargs...))
"""
macro measure(d,b)
    d = esc(d)
    b = esc(b)
    return quote
        struct $d{P,X} <: AbstractMeasure{X}
            par :: P
            $d(nt::NamedTuple) = new{typeof(nt), _domain($d, typeof(nt))}(nt)
        end

        $d(;kwargs...) = $d((;kwargs...))

        baseMeasure($d) = $b

        ≪(::$d{P,X}, ::$b{X}) where {P,X} = true
    end    
end

include("basemeasures.jl")
include("basemeasures/lebesgue.jl")
include("combinators/scale.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("distributions.jl")
# include("probability/normal.jl")



end # module
