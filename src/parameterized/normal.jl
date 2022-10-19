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

massof(::Normal) = static(1.0)

for N in AFFINEPARS
    @eval begin
        proxy(d::Normal{$N,T}) where {T} = affine(params(d), Normal())
        @useproxy Normal{$N,T} where {T}
    end
end

insupport(d::Normal, x) = true

insupport(d::Normal) = Returns(true)

@inline logdensity_def(d::Normal{()}, x) = -x^2 / 2
@inline basemeasure(::Normal{()}) = WeightedMeasure(static(-0.5 * log2π), LebesgueBase())

@kwstruct Normal(μ)
@kwstruct Normal(σ)
@kwstruct Normal(μ, σ)
@kwstruct Normal(λ)
@kwstruct Normal(μ, λ)

params(::Type{N}) where {N<:Normal} = ()

Normal(μ::M, σ::S) where {M,S} = Normal((μ = μ, σ = σ))::Normal{(:μ, :σ),Tuple{M,S}}

# Normal(nt::NamedTuple{N,Tuple{Vararg{AbstractArray}}}) where {N} = MvNormal(nt)

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
    ϕ     => σ²
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

# Rather than try to reimplement everything in Distributions, measures can have
# a `proxy` method. This just delegates some methods to the corresponding
# Distributions.jl methods. For example,
#
#     julia> entropy(Normal(2,4))
#     2.805232894324563
#
proxy(d::Normal{()}) = Dists.Normal()

###############################################################################
# Some distributions have a "standard" version that takes no parameters
@kwstruct Normal()

# Instead of setting default values, the `@kwstruct` call above makes a
# parameter-free instance available. The log-density for this is very efficient.

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::Normal{()}) where {T} = randn(rng, T)
Base.rand(rng::Random.AbstractRNG, ::Type{T}, μ::Normal) where {T} = rand(rng, T, proxy(μ))

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
HalfNormal(σ) = HalfNormal((σ = σ,))

###############################################################################
@kwstruct Normal(σ²)
@kwstruct Normal(μ, σ²)

@inline function logdensity_def(d::Normal{(:σ²,)}, x)
    -0.5 * x^2 / d.σ²
end

@inline function basemeasure(d::Normal{(:σ²,)})
    ℓ = static(-0.5) * (static(float(log2π)) + log(d.σ²))
    weightedmeasure(ℓ, LebesgueBase())
end

proxy(d::Normal{(:μ, :σ²)}) = affine((μ = d.μ,), Normal((σ² = d.σ²,)))
@useproxy Normal{(:μ, :σ²)}

###############################################################################
@kwstruct Normal(τ)
@kwstruct Normal(μ, τ)

@inline function logdensity_def(d::Normal{(:τ,)}, x)
    -0.5 * x^2 * d.τ
end

@inline function basemeasure(d::Normal{(:τ,)})
    ℓ = static(-0.5) * (static(float(log2π)) - log(d.τ))
    weightedmeasure(ℓ, LebesgueBase())
end

proxy(d::Normal{(:μ, :τ)}) = affine((μ = d.μ,), Normal((τ = d.τ,)))
@useproxy Normal{(:μ, :τ)}

###############################################################################
@kwstruct Normal(μ, logσ)

@inline function logdensity_def(d::Normal{(:μ, :logσ)}, x)
    μ = d.μ
    logσ = d.logσ
    - 0.5(exp(-2 * logσ) * ((x - μ)^2))
end

function basemeasure(d::Normal{(:μ, :logσ)})
    ℓ = static(-0.5) * (static(float(log2π)) + static(2.0) * d.logσ)
    weightedmeasure(ℓ, LebesgueBase())
end

function logdensity_def(p::Normal, q::Normal, x)
    μp = mean(p)
    σp = std(p)
    μq = mean(q)
    σq = std(q)

    # Result is (((x - μq) / σq)^2 - ((x - μp) / σp)^2 + log(abs(σq / σp))) / 2 

    # We'll write the difference of squares as sqdiff, then divide that by 2 at
    # the end

    if σp == σq
        return (2x - μq - μp) * (μp - μq) / (2 * σp^2)
    else
        zp = (x - μp) / σp
        zq = (x - μq) / σq
        return ((zq + zp) * (zq - zp)) / 2 + log(abs(σq / σp))
    end
end

MeasureBase.transport_origin(::Normal) = StdNormal()

MeasureBase.to_origin(::Normal{()}, y) = y
MeasureBase.from_origin(::Normal{()}, x) = x

MeasureBase.smf(::Normal{()}, x) = Φ(x)
MeasureBase.invsmf(::Normal{()}, p) = Φinv(p)

@smfAD Normal{()}