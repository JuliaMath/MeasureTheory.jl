module Measures

using MLStyle
import Distributions
using Reexport
using SimplePosets
using SimpleTraits

import Distributions
const Dists = Distributions
# include("traits.jl")

const EmptyNamedTuple = NamedTuple{(),Tuple{}}

abstract type AbstractMeasure{X} end

"""
    logdensity(μ::Measure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `baseMeasure(μ)`.

Methods for computing density relative to other measures will be 
"""
function logdensity end


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
        Normal(nt::NamedTuple) = new{typeof(nt), eltype(Normal, typeof(nt))}(nt)
    end

    Normal(; kwargs...) = Normal((;kwargs...))
"""
macro measure(d,b)
    d = esc(d)
    b = esc(b)
    return quote
        struct $d{P,X} <: AbstractMeasure{X}
            par :: P
            $d(nt::NamedTuple) = new{typeof(nt), eltype($d, typeof(nt))}(nt)
        end

        $d(;kwargs...) = $d((;kwargs...))

        baseMeasure(μ::$d{P,X}) where {P,X} = $b{X}

        ≪(::$d{P,X}, ::$b{X}) where {P,X} = true
    end    
end

include("basemeasures.jl")
include("basemeasures/lebesgue.jl")
include("combinators/scale.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("distributions.jl")
include("probability/normal.jl")



end # module
