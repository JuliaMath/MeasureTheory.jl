module MeasureTheory

using MLStyle

import Distributions
const Dists = Distributions

const EmptyNamedTuple = NamedTuple{(),Tuple{}}


export ≪
export sampletype

abstract type AbstractMeasure end

"""
    logdensity(μ::Measure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `baseMeasure(μ)`.

Methods for computing density relative to other measures will be 
"""
function logdensity end

include("absolutecontinuity.jl")
include("basemeasures.jl")
include("macros.jl")
include("combinators/scale.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("distributions.jl")
include("probability/normal.jl")
include("density.jl")
include("pushforward.jl")

end # module
