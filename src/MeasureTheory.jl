module MeasureTheory

using Compat
using Random

using MeasureBase
using MLStyle

using DistributionMeasures
# using TransformVariables: asℝ₊, as𝕀, asℝ, transform

import Base
import Distributions
const Dists = Distributions

export transform
export ≪
export gentype
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

using Infinities
using DynamicIterators
using KeywordCalls
using ConstructionBase
using Accessors
using StatsFuns
using SpecialFunctions
using ConcreteStructs

import LogExpFunctions
import NamedTupleTools
import InverseFunctions: inverse
export inverse

import MeasureBase: insupport, instance, marginals
import MeasureBase:
    testvalue,
    logdensity_def,
    density_def,
    basemeasure,
    kernel,
    params,
    paramnames,
    ∫,
    𝒹,
    ∫exp
import MeasureBase: ≪
using MeasureBase: BoundedInts, BoundedReals, CountingMeasure, IntegerDomain, IntegerNumbers
using MeasureBase: weightedmeasure, restrict
using MeasureBase: AbstractTransitionKernel

import Statistics: mean, var, std

import MeasureBase: likelihoodof
export likelihoodof
export log_likelihood_ratio

using StaticArraysCore

import PrettyPrinting

const Pretty = PrettyPrinting

import Base: rand

using Reexport
@reexport using MeasureBase
import IfElse
using IfElse

using Tricks: static_hasmethod

using Static

export as
export AffinePushfwd
export AffineTransform
export insupport
export For

using MeasureBase: kernel
using MeasureBase: Returns
import MeasureBase: proxy, @useproxy
import MeasureBase: basemeasure_depth
using MeasureBase: LebesgueMeasure

import DensityInterface: logdensityof
import DensityInterface: densityof
import DensityInterface: DensityKind
using DensityInterface

gentype(μ::AbstractMeasure) = typeof(testvalue(μ))

# gentype(μ::AbstractMeasure) = gentype(basemeasure(μ))

xlogx(x::Number) = LogExpFunctions.xlogx(x)
xlogx(x, y) = x * log(x)

xlogy(x::Number, y::Number) = LogExpFunctions.xlogy(x, y)
xlogy(x, y) = x * log(y)

xlog1py(x::Number, y::Number) = LogExpFunctions.xlog1py(x, y)
xlog1py(x, y) = x * log(1 + y)

as(args...; kwargs...) = TV.as(args...; kwargs...)

include("utils.jl")
# include("const.jl")
include("combinators/for.jl")
# include("traits.jl")
include("parameterized.jl")

include("macros.jl")
include("combinators/affine.jl")
include("combinators/weighted.jl")
include("combinators/product.jl")
# include("combinators/transforms.jl")
# include("combinators/exponential-families.jl")
include("resettable-rng.jl")
include("realized.jl")
include("combinators/chain.jl")

# include("distributions.jl")
include("smart-constructors.jl")

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
include("parameterized/betabinomial.jl")
include("parameterized/gamma.jl")
include("parameterized/snedecorf.jl")
include("parameterized/inverse-gaussian.jl")

include("combinators/ifelse.jl")

# include("transforms/corrcholesky.jl")
# include("transforms/ordered.jl")

include("distproxy.jl")
end # module
