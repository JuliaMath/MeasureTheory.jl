# Adding a New Measure

## Parameterized Measures

This is by far the most common kind of measure, and is especially useful as a way to describe families of proibability distributions.

### Declaring a Parameterized Measure

To start, declare a `@measure`. For example, `Normal` is declared as

```julia
@measure Normal(μ,σ) ≪ (1/sqrt2π) * Lebesgue(ℝ)
```

[`ℝ` is typed as `\bbR <TAB>`]

A `ParameterizedMeasure` can have multiple parameterizations, which as dispatched according to the names of the parameters. The `(μ,σ)` here specifies names to assign if none are given. So for example,

```julia
julia> Normal(-3.0, 2.1)
Normal(μ = -3.0, σ = 2.1)
```

The right side, `(1/sqrt2π) * Lebesgue(ℝ)`, gives the _base measure_. `Lebesgue` in this case is the technical name for the measure associating an interval of real numbers to its length. The `(1/sqrt2π)` comes from the normalization constant in the probability density function,

```math
f_{\text{Normal}}(x|μ,σ) = \frac{1}{σ \sqrt{2 \pi}} e^{-\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2}\ \ .
```

Making this part of the base measure allows us to avoid including it in every computation.

The `≪` (typed as `\ll <TAB>`) can be read "is dominated by". This just means that any set for which the base measure is zero must also have zero measure in what we're defining.

### Defining a Log Density

Most computations involve log-densities rather than densities themselves, so these are our first priority. `density(d,x)` will default to `exp(logdensity(d,x))`, but you can add a separate method if it's more efficient.

The definition is simple:
```julia
logdensity(d::Normal{()} , x) = - x^2 / 2 
```

There are a few things here worth noting.

First, we dispatch by the names of `d` (here there are none), and the type of `x` is not specified.

Also, there's nothing here about `μ` and `σ`. These _location-scale parameters_ behave exactly the same across lots of distributions, so we have a macro to add them:

```julia
@μσ_methods Normal()
```

Let's look at another example, the Beta distribution. Here the base measure is `Lebesgue(𝕀)` (support is the unit interval). The log-density is

```julia
function logdensity(d::Beta{(:α, :β)}, x)
    return (d.α - 1) * log(x) + (d.β - 1) * log(1 - x) - logbeta(d.α, d.β)
end
```

Note that when possible, we avoid extra control flow for checking that `x` is in the support. In applications, log-densities are often evaluated only on the support by design. Such checks should be implemented at a higher level whenever there is any doubt.

Finally, note that in both of these examples, the log-density has a relatively direct algebraic form. It's imnportant to have this whenever possible to allow for symbolic manipulation (using libraries like `SymolicUtils.jl`) and efficient automatic differentiation.

### Random Sampling

For univariate distributions, you should define a `Base.rand` method that uses three arguments:
- A `Random.AbstractRNG` to use for randomization,
- A type to be returned, and
- The measure to sample from.

For our `Normal` example, this is

```julia
Base.rand(rng::Random.AbstractRNG, T::Type, d::Normal{()}) = randn(rng, T)
```

Again, for location-scale families, other methods are derived automatically. 

For multivariate distributions (or anything that requires heap allocation), you should instead define a `Random.rand!` method. This also takes three arguments, different from `rand`:
- The `Random.AbstractRNG`,
- The measure to sample from, and
- Where to store the result.

For example, here's the implementation for `ProductMeasure`:

```julia
@propagate_inbounds function Random.rand!(rng::AbstractRNG, d::ProductMeasure, x::AbstractArray)
    @boundscheck size(d.data) == size(x) || throw(BoundsError)

    @inbounds for j in eachindex(x)
        x[j] = rand(rng, eltype(x), d.data[j])
    end
    x
end
```

Note that in this example, `d.data[j]` might itself require allocation. This implementation is likely to change in the future to make it easier to define nested structures without any need for ongoing allocation.

## Primitive Measures

Most measures are defined in terms of a `logdensity` with respect to some other measure, its `basemeasure`. But how is _that_ measure defined? It can't be "densities all the way down"; at some point, the chain has to stop.

A _primitive_ measure is just a measure that has itself as its own base measure. Note that this also means its log-density is always zero.

Here's the implementation of `Lebesgue`:

```julia
struct Lebesgue{X} <: AbstractMeasure end

Lebesgue(X) = Lebesgue{X}()

basemeasure(μ::Lebesgue) = μ

isprimitive(::Lebesgue) = true

sampletype(::Lebesgue{ℝ}) = Float64
sampletype(::Lebesgue{ℝ₊}) = Float64
sampletype(::Lebesgue{𝕀}) = Float64

logdensity(::Lebesgue, x) = zero(float(x))
```
