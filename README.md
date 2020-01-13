# Measures

[![Build Status](https://travis-ci.com/cscherrer/Measures.jl.svg?branch=master)](https://travis-ci.com/cscherrer/Measures.jl)
[![Codecov](https://codecov.io/gh/cscherrer/Measures.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cscherrer/Measures.jl)

`Measures.jl` is a package for building and reasoning about measures.

For an example, let's walk through the construction of `src/probability/Normal`.

First, we have

```julia
@measure Normal
```

this is just a little helper function, and is equivalent to

```julia
struct Normal{P, X} <: Measure{X}
    par::P
    Normal(nt::NamedTuple) = new{typeof(nt), _domain(Normal, typeof(nt))}(nt)
end

Normal(; kwargs...) = Normal((;kwargs...))
```

Next we have 

```julia
Normal(μ::Real, σ::Real) = Normal(μ=μ, σ=σ)
```

This defines a default. If we just give two numbers as arguments (but no names), we'll assume this parameterization.

Next need to define a `_domain` function. This takes a constructor (here `Normal`) and a parameter, and tells us the space for which this defines a measure. Let's define this in terms of the types of the parameters,

```julia
_domain(::Type{Normal}, ::Type{NamedTuple{(:μ, :σ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)
```

That's still kind of boring, so let's build the density. For this, we need to implement the trait

```julia
@trait Density{M,X} where {X = domain{M}} begin
    baseMeasure :: [M] => Measure{X}
    logdensity :: [M, X] => Real
end
```

A density doesn't exist by itself, but is defined relative to some _base measure_. For a normal distribution this is just Lebesgue measure on the real numbers. That, together with the usual Gaussian log-density, gives us

```julia
@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :σ)}} begin
    baseMeasure(d) = Lebesgue(X)
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
_domain(::Type{Normal}, ::Type{NamedTuple{(:μ, :τ), Tuple{A, B}}}) where {A,B} = promote_type(A,B)

@implement Density{Normal{P,X},X} where {X, P <: NamedTuple{(:μ, :τ)}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d, x) = - (log(2) + log(π) - log(d.par.τ)  + d.par.τ * (x - d.par.μ)^2) / 2
end
```

And another check:

```julia
julia> logdensity(Normal(μ=0.0, τ=4.0), 1.0)
-2.2257913526447273
```
