# Variate Transformations

Often one needs to find a map `f` that transports samples of probability measure `μ` into samples of another probability measure `v`. More generally, `μ` to `v` can be measures. Then `v` is the [pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure) of  measure `μ` under such `f`.

Using [`vartransform`](@ref) we can (for many measures) automatically generate this function `f = vartransform(v, μ)` so that [`pushfwd(f, μ)`](@ref) becomes equal to `v`. Instead of sampling `∫(L, v)` we can now sample `∫(L∘f, μ)` by transforming the sampled points `X_μ`  to sample points `X_v = f.(X_μ)`.

A typical application arises in sampling: Many sampling algorithms perform better or even require a specific  base measure. Nested sampling, for example, natively uses μ as it's base measure, while many MCMC sampling algorithms work best in an unbounded space and prefer base measures such as `μ = StdNormal()^n`. Using [`vartransform`](@ref) we can (for many measures) automatically generate a function `f = vartransform(v, μ)` so that [`pushfwd(f, μ)`](@ref) becomes equivalent to `v`. Instead of sampling `∫(L, v)` we can now sample `∫(L∘f, μ)`. The generates sample points `X_μ` can then be transformed to sample points `X_v = f.(X_μ)`.

Transformation functions `f = vartransform(v, μ)` support `f_inverse = InverseFunctions.inverse(f)` and `x_v, ladj = ChangesOfVariables.with_logabsdet_jacobian(f, x_μ)`.
