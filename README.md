# MeasureTheory

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cscherrer.github.io/MeasureTheory.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cscherrer.github.io/MeasureTheory.jl/dev)
[![Build Status](https://github.com/cscherrer/MeasureTheory.jl/workflows/CI/badge.svg)](https://github.com/cscherrer/MeasureTheory.jl/actions)
[![Coverage](https://codecov.io/gh/cscherrer/MeasureTheory.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cscherrer/MeasureTheory.jl)

`MeasureTheory.jl` is a package for building and reasoning about measures.

# Why?

A distribution (as in Distributions.jl) is also called a _probability measure_, and carries with it the constraint of adding (or integrating) to one. Statistical work usually requires this "at the end of the day", but enforcing it at each step of a computation can have considerable overhead.

As a generalization of the concept of volume, measures also have applications outside of probability theory.

# Goals

## Distributions.jl Compatibility

Distirbutions.jl is wildly popular, and is large enough that replacing it all at once would be a major undertaking. 

Instead, we should aim to make any Distribution easily usable as a Measure. We'll most likely implement this using an `IsMeasure` trait. 

## Absolute Continuity

For two measures μ, ν on a set X, we say μ is _absolutely continuous_ with respect to ν if ν(A)=0 implies μ(A)=0 for every measurable subset A of X.

The following are equivalent:
1. μ ≪ ν
2. μ is absolutely continuous wrt ν
3. There exists a function f such that μ = ∫f dν

So we'll need a `≪` operator. Note that `≪` is not antisymmetric; it's common for both `μ ≪ ν` and  `ν ≪ μ` to hold. 

If `μ ≪ ν` and  `ν ≪ μ`, we say μ and ν are _equivalent_ and write `μ ≃ ν`. (This is often written as `μ ~ ν`, but we reserve `~` for random variables following a distribution, as is common in the literature and probabilistic programming languages.)

If we collapse the equivalence classes (under ≃), `≪` becomes a partial order.

_We need ≃ and ≪ to be fast_. If the support of a measure can be determined statically from its type, we can define ≃ and ≪ as generated functions. 

## Radon-Nikodym Derivatives

One of the equivalent conditions above was "There exists a function f such that μ = ∫f dν". In this case, `f` is called a _Radon-Nikodym derivative_, or (less formally) a _density_. In this case we often write `f = dμ/dν`.

For any measures μ and ν with μ≪ν, we should be able to represent this.

## Integration

More generally, we'll need to be able to represent change of measure as above, `∫f dν`. We'll need an `Integral` type

```julia
struct Integral{F,M}
    f::F
    μ::M
end
```

Then we'll have a function `∫`. For cases where μ = ∫f dν,  `∫(f, ν)` will just return `μ` (we can do this based on the types). For unknown cases (which will be most of them), we'll return `∫(f, ν) = Integral(f, ν)`.

## Measure Combinators

It should be very easy to build new measures from existing ones. This can be done using, for example, 

- restriction
- product measure
- superposition
- pushforward

There's also function spaces, but this gets much trickier. We'll need to determine a good way to reason about this.

## More???

This is very much a work in progress. If there are things you think we should have as goals, please add an issue with the `goals` label.


------------------
# Old Stuff

**WARNING: The current README is very developer-oriented. Casual use will be much simpler**

For an example, let's walk through the construction of `src/probability/Normal`.

First, we have

```julia
@measure Normal
```

this is just a little helper function, and is equivalent to

# TODO: Clean up
```julia
quote
    #= /home/chad/git/Measures.jl/src/Measures.jl:55 =#
    struct Normal{var"#10#P", var"#11#X"} <: Measures.AbstractMeasure{var"#11#X"}
        #= /home/chad/git/Measures.jl/src/Measures.jl:56 =#
        par::var"#10#P"
    end
    #= /home/chad/git/Measures.jl/src/Measures.jl:59 =#
    function Normal(var"#13#nt"::Measures.NamedTuple)
        #= /home/chad/git/Measures.jl/src/Measures.jl:59 =#
        #= /home/chad/git/Measures.jl/src/Measures.jl:60 =#
        var"#12#P" = Measures.typeof(var"#13#nt")
        #= /home/chad/git/Measures.jl/src/Measures.jl:61 =#
        return Normal{var"#12#P", Measures.eltype(Normal{var"#12#P"})}
    end
    #= /home/chad/git/Measures.jl/src/Measures.jl:64 =#
    Normal(; var"#14#kwargs"...) = begin
            #= /home/chad/git/Measures.jl/src/Measures.jl:64 =#
            Normal((; var"#14#kwargs"...))
        end
    #= /home/chad/git/Measures.jl/src/Measures.jl:66 =#
    (var"#8#basemeasure"(var"#15#μ"::Normal{var"#16#P", var"#17#X"}) where {var"#16#P", var"#17#X"}) = begin
            #= /home/chad/git/Measures.jl/src/Measures.jl:66 =#
            Lebesgue{var"#17#X"}
        end
    #= /home/chad/git/Measures.jl/src/Measures.jl:68 =#
    (var"#9#≪"(::Normal{var"#19#P", var"#20#X"}, ::Lebesgue{var"#20#X"}) where {var"#19#P", var"#20#X"}) = begin
            #= /home/chad/git/Measures.jl/src/Measures.jl:68 =#
            true
        end
end
```

Next we have 

```julia
Normal(μ::Real, σ::Real) = Normal(μ=μ, σ=σ)
```

This defines a default. If we just give two numbers as arguments (but no names), we'll assume this parameterization.

Next need to define a `eltype` function. This takes a constructor (here `Normal`) and a parameter, and tells us the space for which this defines a measure. Let's define this in terms of the types of the parameters,

```julia
eltype(::Type{Normal}, ::Type{NamedTuple{(:μ, :σ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)
```

That's still kind of boring, so let's build the density. For this, we need to implement the trait

```julia
@trait Density{M,X} where {X = domain{M}} begin
    basemeasure :: [M] => Measure{X}
    logdensity :: [M, X] => Real
end
```

A density doesn't exist by itself, but is defined relative to some _base measure_. For a normal distribution this is just Lebesgue measure on the real numbers. That, together with the usual Gaussian log-density, gives us

```julia
@implement Density{Normal{X,P},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
    basemeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π)) / 2 - log(d.par.σ)  - (x - d.par.μ)^2 / (2 * d.par.σ^2)
end
```

Now we can compute the log-density:

```julia
julia> logdensity(Normal(0.0, 0.5), 1.0)
-2.2257913526447273
```

And just to check that our default is working,

```julia
julia> logdensity(Normal(μ=0.0, σ=0.5), 1.0)
-2.2257913526447273
```

What about other parameterizations? Sure, no problem. Here's a way to write this for mean `μ` (as before), but using the _precision_ (inverse of the variance) instead of standard deviation:

```julia
eltype(::Type{Normal}, ::Type{NamedTuple{(:μ, :τ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)

@implement Density{Normal{X,P},X} where {X, P <: NamedTuple{(:μ, :τ)}} begin
    basemeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π) - log(d.par.τ)  + d.par.τ * (x - d.par.μ)^2) / 2
end
```

And another check:

```julia
julia> logdensity(Normal(μ=0.0, τ=4.0), 1.0)
-2.2257913526447273
```

We can combine measures in a few ways, for now just _scaling_ and _superposition_:

```julia
julia> 2.0*Lebesgue(Float64) + Normal(0.0,1.0)
SuperpositionMeasure{Float64,2}((MeasureTheory.WeightedMeasure{Float64,Float64}(2.0, Lebesgue{Float64}()), Normal{NamedTuple{(:μ, :σ),Tuple{Float64,Float64}},Float64}((μ = 0.0, σ = 1.0))))
```

---

For an easy way to find expressions for the common log-densities, see [this gist](https://gist.github.com/cscherrer/47f0fc7597b4ffc11186d54cc4d8e577)

## Support

[![Planting Space](https://planting.space/sponsor/PlantingSpace-sponsor-3.png)](https://planting.space)

## Stargazers over time

[![Stargazers over time](https://starchart.cc/cscherrer/MeasureTheory.jl.svg)](https://starchart.cc/cscherrer/MeasureTheory.jl)
