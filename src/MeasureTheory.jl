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

export TV
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
export basekernel

using InfiniteArrays
using ConcreteStructs
using DynamicIterators
using KeywordCalls
using ConstructionBase
using Accessors
using StatsFuns
using SpecialFunctions

import LogExpFunctions 
import NamedTupleTools

import MeasureBase: testvalue, logdensity, density, basemeasure, kernel, params, ∫
import MeasureBase: affine, supportdim

import Base: rand

using Reexport
@reexport using MeasureBase

using Tricks: static_hasmethod
const ∞ = InfiniteArrays.∞

export ∞

export as
export Affine
export AffineTransform

using MeasureBase: Returns

sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

import Distributions: logpdf, pdf

export pdf, logpdf

Distributions.logpdf(d::AbstractMeasure, x) = MeasureBase.logpdf(d, x)

Distributions.pdf(d::AbstractMeasure, x) = exp(Dists.logpdf(d, x))

"""
    logdensity(μ::AbstractMeasure [, ν::AbstractMeasure], x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.
"""
function logdensity end


const AFFINEPARS = [
    (:μ,:σ)
    (:μ,:ω)
    (:σ,)
    (:ω,)
    (:μ,)
]

xlogy(x::Number, y::Number) = LogExpFunctions.xlogy(x, y)
xlogy(x, y) = x * log(y)

xlog1py(x::Number, y::Number) = LogExpFunctions.xlog1py(x, y)
xlog1py(x, y) = x * log1p(y)


include("const.jl")
# include("traits.jl")
include("parameterized.jl")
# include("resettablerng.jl")

include("combinators/affine.jl")
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
