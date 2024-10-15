module MeasureTheory

using Compat
using Random

using MeasureBase
using MLStyle
import TransformVariables
const TV = TransformVariables

using TransformVariables: asℝ₊, as𝕀, asℝ, transform

import Base

export TV
export transform
export For

export AbstractMeasure
export Dirac
export Lebesgue
export ℝ, ℝ₊, 𝕀
export ⊙
export SpikeMixture
export TrivialMeasure
export Likelihood
export testvalue
export basekernel

using Infinities
using KeywordCalls
using ConstructionBase
using Accessors
# using StatsFuns
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
    ∫exp,
    smf,
    invsmf,
    massof
using MeasureBase: BoundedInts, BoundedReals, CountingBase, IntegerDomain, IntegerNumbers
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
using MeasureBase: LebesgueBase
using ConstantRNGs

import DensityInterface: logdensityof
import DensityInterface: densityof
import DensityInterface: DensityKind
using DensityInterface

using ForwardDiff
using ForwardDiff: Dual


xlogx(x::AbstractFloat) = LogExpFunctions.xlogx(x)
xlogx(x, y) = x * log(x)

xlogy(x::AbstractFloat, y::AbstractFloat) = LogExpFunctions.xlogy(x, y)
xlogy(x, y) = x * log(y)

xlog1py(x::AbstractFloat, y::AbstractFloat) = LogExpFunctions.xlog1py(x, y)
xlog1py(x, y) = x * log(1 + y)

log1pexp(x::AbstractFloat) = LogExpFunctions.log1pexp(x)
log1pexp(x) = log(1 + exp(x))

using MeasureBase: Φ, Φinv
as(args...; kwargs...) = TV.as(args...; kwargs...)

include("consts.jl")
include("utils.jl")
include("const.jl")
include("combinators/for.jl")

include("macros.jl")
include("combinators/affine.jl")
include("combinators/weighted.jl")
include("combinators/transforms.jl")
# include("combinators/exponential-families.jl")
# TODO: Clean these up and re-add them
# include("resettable-rng.jl")
# include("realized.jl")
# include("combinators/chain.jl")

include("distributions.jl")
include("smart-constructors.jl")

include("scalefree.jl")

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
include("parameterized/inverse-gamma.jl")
include("parameterized/bernoulli.jl")
include("parameterized/poisson.jl")
include("parameterized/binomial.jl")
include("parameterized/multinomial.jl")
include("parameterized/lkj-cholesky.jl")
include("parameterized/negativebinomial.jl")
include("parameterized/betabinomial.jl")
include("parameterized/gamma.jl")
include("parameterized/snedecorf.jl")
# include("parameterized/inverse-gaussian.jl")
include("parameterized/scaled-invchisq.jl")
include("parameterized/normal-invchisq.jl")

include("combinators/ifelse.jl")

include("transforms/corrcholesky.jl")
include("transforms/ordered.jl")

include("parameterized.jl")

end # module
