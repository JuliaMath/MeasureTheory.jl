module MeasureTheory

using Random

using MeasureBase
using ConcreteStructs
using MLStyle
using NestedTuples
using TransformVariables
const TV = TransformVariables

import Base
import Distributions
const Dists = Distributions

export ≪
export sampletype
export For

export AbstractMeasure
export Dirac
export Lebesgue
export ℝ, ℝ₊, 𝕀
export ⊙
export SpikeMixture
export CountingMeasure
export TrivialMeasure
export Likelihood
export testvalue

using InfiniteArrays
using ConcreteStructs
using DynamicIterators
using KeywordCalls
using ConstructionBase
using Accessors
using StatsFuns
using SpecialFunctions

import MeasureBase: testvalue, logdensity, density, basemeasure, kernel
using Tricks: static_hasmethod
const ∞ = InfiniteArrays.∞

export ∞
export testvalue, logdensity, density, basemeasure, kernel

export as
export Affine
export AffineTransform

sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

"""
    logdensity(μ::AbstractMeasure [, ν::AbstractMeasure], x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.
"""
function logdensity end

include("const.jl")
# include("traits.jl")
include("parameterized.jl")
# include("resettablerng.jl")

include("combinators/weighted.jl")
include("combinators/product.jl")
include("combinators/transforms.jl")
include("combinators/chain.jl")

include("distributions.jl")

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
include("parameterized/lkj-cholesky.jl")
include("parameterized/negativebinomial.jl")

include("transforms/corrcholesky.jl")
include("transforms/ordered.jl")

include("distproxy.jl")
end # module
