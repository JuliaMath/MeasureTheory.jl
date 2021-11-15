using StatsFuns
export Normal, HalfNormal

# The `@parameterized` macro below does three things:
# 1. Defines the `Normal` struct 
#
#    struct Normal{N, T} <: (ParameterizedMeasure){N}
#        par::MeasureTheory.NamedTuple{N, T}
#    end
#
#    (Note that this is the same form required by the `@kwstruct` in the
#    KeywordCalls package.)  
#
# 2. Sets the base measure
#
#    basemeasure(::Normal) = (1 / sqrt2π) * Lebesgue(ℝ)
#
#    This just means that (log-)densities will be defined relative to this.
#    Including the `(1 / sqrt2π)` in the base measure allows us to avoid
#    carrying around the normalization constant.
#   
# 3. Sets up a method for a call with two unnamed arguments,
#
#    Normal(μ, σ) = Normal((μ=μ, σ=σ))
#
@parameterized Normal() 

basemeasure(::Normal{()}) = WeightedMeasure(-0.5*log2π, Lebesgue(ℝ))

basemeasure_depth(::Type{Normal{()}}) = static(2)

basemeasure_depth(::Type{<:Normal}) = static(3)

@kwstruct Normal(μ)
@kwstruct Normal(σ)
@kwstruct Normal(μ,σ)
@kwstruct Normal(ω)
@kwstruct Normal(μ,ω)

params(::Type{N}) where {N<:Normal} = ()

Normal(μ,σ) = Normal((μ=μ, σ=σ))

using MeasureBase: rowsize, colsize

Normal(nt::NamedTuple{N,Tuple{Vararg{AbstractArray}}}) where {N} = MvNormal(nt)

for N in AFFINEPARS
    @eval begin
        proxy(d::Normal{$N}) = affine(params(d), Normal())
        logdensity(d::Normal{$N}, x) = logdensity(proxy(d), x)
        basemeasure(d::Normal{$N}) = basemeasure(proxy(d))
    end
end

TV.as(::Normal) = asℝ

# `@kwalias` defines some alias names, giving users flexibility in the names
# they use. For example, σ² is standard notation for the variance parameter, but
# it's a lot to type. Some users might prefer to just use `var` and have us do
# the conversion (at compile time).
@kwalias Normal [
    mean  => μ
    mu    => μ
    std   => σ
    sigma => σ
    var   => σ²
]

# It's often useful to be able to map into the parameter space for a given
# measure. We (currently, at least) do this *independently* per parameter. This
# allows us to do things like, e.g.
#
# julia> transform(asparams(Normal(2,3)))(randn(2))
# (μ = -0.6778046205440658, σ = 3.7357686957861898)
#
# julia> transform(asparams(Normal(μ=-3,logσ=2)))(randn(2))
# (μ = 0.3995295982002209, logσ = -1.3902312393777492)
# 
# Or using types:
# julia> transform(asparams(Normal{(:μ,:σ²)}))(randn(2))
# (μ = -0.4548087051528626, σ² = 11.920775478312793)
#
# And of course, you can apply `Normal` to any one of the above.
#
asparams(::Type{<:Normal}, ::Val{:σ²}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:τ}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:logτ}) = asℝ


# Rather than try to reimplement everything in Distributions, measures can have
# a `distproxy` method. This just delegates some methods to the corresponding
# Distributions.jl methods. For example,
#
#     julia> entropy(Normal(2,4))
#     2.805232894324563
#
distproxy(d::Normal{()}) = Dists.Normal()
distproxy(d::Normal{(:μ,)}) = Dists.Normal(d.μ, 1.0)
distproxy(d::Normal{(:σ,)}) = Dists.Normal(0.0, d.σ)
distproxy(d::Normal{(:μ,:σ)}) = Dists.Normal(d.μ, d.σ)
distproxy(d::Normal{(:ω,)}) = Dists.Normal(0.0, inv(d.ω))
distproxy(d::Normal{(:μ,:ω)}) = Dists.Normal(d.μ, inv(d.ω))


###############################################################################
# Some distributions have a "standard" version that takes no parameters
@kwstruct Normal()

# Instead of setting default values, the `@kwstruct` call above makes a
# parameter-free instance available. The log-density for this is very efficient.
logdensity(d::Normal{()} , x) = - x^2 / 2 

Base.rand(rng::Random.AbstractRNG, T::Type, μ::Normal{()}) = randn(rng, T)

###############################################################################
# `μ` and `σ` parameterizations are so common, we have a macro to make them easy
# to build. Note that `μ` and `σ` are *not* always mean and standard deviation,
# but instead should be considered "location" and "scale" parameters,
# respectively. Also, it's valid to have existing parameters here as a starting
# point, in particular a "shape parameter" is very common. See `StudentT` for an
# example of this.
#
# This defines methods for `logdensity` and `Base.rand`. The latter works out
# especially nicely, because it lets us do things like
#    
#     julia> using Symbolics
#     [ Info: Precompiling Symbolics [0c5d862f-8b57-4792-8d23-62f2024744c7]
#    
#     julia> @variables μ σ
#     2-element Vector{Num}:
#      μ
#      σ
#    
#     julia> rand(Normal(μ,σ))
#     μ + 1.2517620265570832σ
#
# @μσ_methods Normal()

###############################################################################
# The `@half` macro takes a symmetric univariate measure and efficiently creates
# a truncated version. 
@half Normal


# A single unnamed parameter for `HalfNormal` should be interpreted as a `σ`
HalfNormal(σ) = HalfNormal(σ = σ)


###############################################################################
@kwstruct Normal(μ,σ²)

@inline function logdensity(d::Normal{(:σ²)}, x)
    σ² = d.σ²
    -0.5 * (log(σ²) + (x^2/σ²))
end

@inline function logdensity(d::Normal{(:μ,:σ²)}, x)
    μ = d.μ
    σ² = d.σ²
    -0.5 * (log(σ²) + ((x - μ)^2/σ²))
end

###############################################################################
@kwstruct Normal(μ,τ)

@inline function logdensity(d::Normal{(:τ)}, x)
    τ = d.τ
    0.5 * (log(τ) - τ * x^2)
end

@inline function logdensity(d::Normal{(:μ,:τ)}, x)
    μ = d.μ
    τ = d.τ
    0.5 * (log(τ) - τ * (x - μ)^2)
end


###############################################################################
@kwstruct Normal(μ, logσ)

@inline function logdensity(d::Normal{(:μ,:logσ)}, x)
    μ = d.μ
    logσ = d.logσ
    -logσ - 0.5(exp(-2logσ)*((x - μ)^2))
end
