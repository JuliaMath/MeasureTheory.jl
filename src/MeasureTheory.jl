module MeasureTheory

using Random

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
is understood to be `basemeasure(μ)`.

Methods for computing density relative to other measures will be
"""
function logdensity end

include("absolutecontinuity.jl")
include("basemeasures.jl")
include("parameterized.jl")
include("macros.jl")
include("combinators/scale.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/power.jl")
include("distributions.jl")
include("rand.jl")
include("probability/dirac.jl")
include("probability/normal.jl")
include("density.jl")
include("pushforward.jl")
include("kernel.jl")
end # module
