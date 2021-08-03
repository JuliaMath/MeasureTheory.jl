module MeasureTheory

using Random

using ConcreteStructs
using MLStyle
using NestedTuples
# using TransformVariables
# const TV = TransformVariables

import Base
import Distributions
const Dists = Distributions

export ≪
export sampletype

export AbstractMeasure
using InfiniteArrays
using ConcreteStructs
using DynamicIterators
using KeywordCalls
using ConstructionBase
using Accessors

const ∞ = InfiniteArrays.∞

export ∞

export as

abstract type AbstractMeasure end

sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

"""
    logdensity(μ::AbstractMeasure [, ν::AbstractMeasure], x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.
"""
function logdensity end

# include("const.jl")
include("exp.jl")
include("domains.jl")
include("utils.jl")
include("traits.jl")
include("absolutecontinuity.jl")
include("parameterized.jl")
include("macros.jl")
include("resettablerng.jl")

include("primitive.jl")
include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/for.jl")
include("combinators/power.jl")
# include("combinators/transforms.jl")
include("combinators/spikemixture.jl")
include("combinators/chain.jl")
include("kernel.jl")
include("combinators/likelihood.jl")
include("combinators/pointwise.jl")

include("distributions.jl")
include("rand.jl")

include("parameterized/normal.jl")
include("parameterized/studentt.jl")
include("parameterized/cauchy.jl")
include("parameterized/laplace.jl")
include("parameterized/uniform.jl")
include("parameterized/beta.jl")
include("parameterized/dirichlet.jl")
include("parameterized/gumbel.jl")
include("parameterized/exponential.jl")
include("parameterized/mvnormal.jl")
# include("parameterized/inverse-gamma.jl")
include("parameterized/bernoulli.jl")
include("parameterized/poisson.jl")
include("parameterized/binomial.jl")
include("parameterized/multinomial.jl")
# include("parameterized/lkj-cholesky.jl")
include("parameterized/negativebinomial.jl")

# include("transforms/corrcholesky.jl")
# include("transforms/ordered.jl")

include("density.jl")
# include("pushforward.jl")

include("distproxy.jl")
end # module
