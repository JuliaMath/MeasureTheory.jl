module MeasureTheory

using Random

using MLStyle
using NestedTuples
using TransformVariables
import Base
import Distributions
const Dists = Distributions

export ≪
export sampletype

export AbstractMeasure
using InfiniteArrays

const ∞ = InfiniteArrays.∞

export ∞

abstract type AbstractMeasure end

sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

"""
    logdensity(μ::AbstractMeasure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.

Methods for computing density relative to other measures will be
"""
function logdensity end

include("paramorder.jl")
include("exp.jl")
include("domains.jl")
include("utils.jl")
include("absolutecontinuity.jl")
include("basemeasures.jl")
include("parameterized.jl")
include("macros.jl")
include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/for.jl")
include("combinators/power.jl")
include("combinators/elementwise.jl")
include("combinators/transforms.jl")
include("combinators/spikemixture.jl")
include("distributions.jl")
include("rand.jl")
include("probability/dirac.jl")
include("probability/normal.jl")
include("probability/studentt.jl")
include("probability/cauchy.jl")
include("probability/laplace.jl")
include("probability/uniform.jl")
include("probability/beta.jl")
include("probability/gumbel.jl")
include("probability/exponential.jl")
include("probability/mvnormal.jl")
include("probability/inverse-gamma.jl")
include("probability/bernoulli.jl")
include("probability/poisson.jl")
include("probability/binomial.jl")
include("probability/LKJL.jl")
include("density.jl")
include("likelihood.jl")
# include("pushforward.jl")
include("kernel.jl")
include("distproxy.jl")
end # module
