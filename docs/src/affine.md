# Affine Transformations

It's very common for measures to use parameters `μ` and `σ`, for example as in `Normal(μ=3, σ=4)` or `StudentT(ν=1, μ=3, σ=4)`. In this context, `μ` and `σ` need not always refer to the mean and standard deviation (the `StudentT` measure specified above is equivalent to a [Cauchy](https://en.wikipedia.org/wiki/Cauchy_distribution) measure, so both mean and standard deviation are undefined).

In general, `μ` is a "location parameter", and `σ` is a "scale parameter". Together these parameters determine an affine transformation.

```math
f(z) = σ z + μ
```

Starting with the above definition, we'll use ``z`` to represent an "un-transformed" variable, typically coming from a measure which has neither a location nor a scale parameter, for example `Normal()`.

It is a common misconception that affine transformations are "linear". Linearity requires ``f(ax + by) = a f(x) + b f(y)`` for scalars ``a`` and ``b``, which would only holds for the above affine transformation ``f`` if ``μ=0``.


## Cholesky-based parameterizations

If the "un-transformed" `z` is univariate, things are relatively simple. But it's important our approach handle the multivariate case as well.

In the literature, it's common for a multivariate normal distribution to be parameterized by a mean `μ` and covariance matrix `Σ`. This is mathematically convenient, but less idea for efficient computation.

While MeasureTheory.jl includes (or will include) a parameterization using `Σ`, we prefer to work in terms of its Cholesky decomposition ``σ``.

The relationship between the computationally efficient "``σ``" and more familiar parameteriation `Σ` can be seen as follows:

Let ``σ`` be a lower-triangular matrix satisfying

```math
σ σᵗ = Σ
```

Then given a (multivariate) standard normal ``z``, the covariance matrix of ``σ z + μ`` is

```math
𝕍[σ z + μ] = Σ
```

The one-dimensional case where we have

```math
𝕍[σ z + μ] = σ²
```

shows that the lower Cholesky factor of the covariance generalizes the concept of standard deviation, completing the link between ``σ`` and `Σ`.

## The "Cholesky precision" parameterization

The ``(μ,σ)`` parameterization is especially convenient for random sampling. Any measure `z ~ Normal()` determines an `x ~ Normal(μ,σ)` through the affine transformation

```math
x = σ z + μ
```

The log-density transformation of a `Normal` with parameters μ, σ does not follow as directly. Starting with an ``x``, we need to find ``z`` using

```math
z = σ⁻¹ (x - μ)
```

so the log-density is

```julia
logdensity(d::Normal{(:μ,:σ)}, x) = logdensity(d.σ \ (x - d.μ)) - logdet(d.σ)
```

Here the `- logdet(σ)` is the "log absolute Jacobian", required to account for the stretching of the space.

The above requires solving a linear system, which adds some overhead. Even with the convenience of a lower triangular system, it's still not quite as efficient as multiplication.

In addition to the covariance ``Σ``, it's also common to parameterize a multivariate normal by its _precision matrix_, defined as the inverse of the covariance matrix, ``Ω = Σ⁻¹``. Similar to our use of ``σ`` for the lower Cholesky factor of `Σ`, we'll use ``ω`` for the lower Cholesky factor of ``Ω``.

This parameterization enables more efficient calculation of the log-density using only multiplication and addition,

```julia
logdensity(d::Normal{(:μ,:ω)}, x) = logdensity(d.ω * (x - d.μ)) + logdet(d.ω)
```

## `AffineTransform`

Transforms like ``z → σ z + μ`` and ``z → ω \ z + μ`` can be specified in MeasureTheory.jl using an `AffineTransform`. For example,

```julia
julia> f = AffineTransform((μ=3.,σ=2.))
AffineTransform{(:μ, :σ), Tuple{Float64, Float64}}((μ = 3.0, σ = 2.0))

julia> f(1.0)
5.0
```

In the univariate case this is relatively simple to invert. But if `σ` is a matrix, matrix inversion becomes necessary. This is not always possible as lower triangular matrices are not closed under matrix inversion and as such are not guaranteed to exist. 

With multiple parameterizations of a given family of measures, we can work around these issues. The inverse transform of a ``(μ,σ)`` transform will be in terms of ``(μ,ω)``, and vice-versa. So

```julia
julia> f⁻¹ = inv(f)
AffineTransform{(:μ, :ω), Tuple{Float64, Float64}}((μ = -1.5, ω = 2.0))

julia> f(f⁻¹(4))
4.0

julia> f⁻¹(f(4))
4.0
```

## `Affine`

Of particular interest (the whole point of all of this, really) is to have a natural way to work with affine transformations of measures. In accordance with the principle of "common things should have shorter names", we call this `Affine`.

The structure of `Affine` is relatively simple:

```julia
struct Affine{N,M,T} <: AbstractMeasure
    f::AffineTransform{N,T}
    parent::M
end
```