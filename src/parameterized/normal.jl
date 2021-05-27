using StatsFuns
export Normal

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
@parameterized Normal(μ,σ) ≪ (1/sqrt2π) * Lebesgue(ℝ)

@kwstruct Normal(μ,Σ)
@kwstruct Normal(μ,Σ⁻¹)

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
    cov   => Σ
]

distproxy(d::Normal{(:μ, :σ)}) = Dists.Normal(d.μ, d.σ)

asparams(::Type{<:Normal}, ::Val{:μ}) = asℝ
asparams(::Type{<:Normal}, ::Val{:σ}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:logσ}) = asℝ
asparams(::Type{<:Normal}, ::Val{:σ²}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:τ}) = asℝ₊
asparams(::Type{<:Normal}, ::Val{:logτ}) = asℝ




Base.rand(rng::Random.AbstractRNG, T::Type, μ::Normal{()}) = randn(rng, T)

###############################################################################
# Some distributions have a "standard" version that takes no parameters
@kwstruct Normal()

# Instead of setting default values, the `@kwstruct` call above makes a
# parameter-free instance available. The log-density for this is very efficient.
logdensity(d::Normal{()} , x) = - x^2 / 2 

###############################################################################
# We need a `@kwstruct` for each instance, including the one declared in the
# `@parameterized` declaration. 
@kwstruct Normal(μ,σ)

# `μ` and `σ` parameterizations are so common, we have a macro to make them easy
# to build. Note that `μ` and `σ` are *not* always mean and standard deviation,
# but instead should be considered "location" and "scale" parameters,
# respectively. Also, it's valid to have existing parameters here as a starting
# point, in particular a "shape parameter" is very common. See `StudentT` for an
# example of this.
@μσ_methods Normal()

###############################################################################
# The `@half` macro takes a symmetric univariate measure and efficiently creates
# a truncated version. 
@half Normal()

@kwstruct HalfNormal()
@kwstruct HalfNormal(σ)

HalfNormal(σ) = HalfNormal(σ = σ)


###############################################################################
@kwstruct Normal(μ,σ²)

function logdensity(d::Normal{(:σ²)}, x)
    σ² = d.σ²
    -0.5 * (log(σ²) + (x^2/σ²))
end

function logdensity(d::Normal{(:μ,:σ²)}, x)
    μ = d.μ
    σ² = d.σ²
    -0.5 * (log(σ²) + ((x - μ)^2/σ²))
end

###############################################################################
@kwstruct Normal(μ,τ)

function logdensity(d::Normal{(:τ)}, x)
    τ = d.τ
    0.5 * (log(τ) - τ * x^2)
end

function logdensity(d::Normal{(:μ,:τ)}, x)
    μ = d.μ
    τ = d.τ
    0.5 * (log(τ) - τ * (x - μ)^2)
end


###############################################################################
@kwstruct Normal(μ, logσ)

function logdensity(d::Normal{(:μ,:logσ)}, x)
    μ = d.μ
    logσ = d.logσ
    -logσ - 0.5(exp(-2logσ)*((x - μ)^2))
end
