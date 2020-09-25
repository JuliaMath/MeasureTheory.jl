module MeasureTheory

using MLStyle

import Distributions
const Dists = Distributions

const EmptyNamedTuple = NamedTuple{(),Tuple{}}


export ≪

abstract type AbstractMeasure{X} end

"""
    logdensity(μ::Measure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `baseMeasure(μ)`.

Methods for computing density relative to other measures will be 
"""
function logdensity end


Base.eltype(μ::AbstractMeasure{X}) where {X} = X



logdensity(μ::Dists.Distribution, x) = Dists.logpdf(μ,x)

density(μ::Dists.Distribution, x) = Dists.pdf(μ,x)


include("absolutecontinuity.jl")
include("basemeasures/lebesgue.jl")
include("macros.jl")
include("combinators/scale.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("distributions.jl")
# include("probability/normal.jl")
include("density.jl")



end # module
