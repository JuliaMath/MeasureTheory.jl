# Variate Transformations

Given a measure `v` that can be seen as a [pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure) of a measure `μ` one often needs to find the [pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure) from `μ` to `v`.

A typical application arises in sampling: Many sampling algorithms perform better or even require a specific  base measure. Nested sampling, for example, natively uses μ as it's base measure, while many MCMC sampling algorithms work best in an unbounded space and prefer base measures such as `μ = StdNormal()^n`. Using [`vartransform`](@ref) we can (for many measures) automatically generate a function `f = vartransform(v, μ)` so that [`pushfwd(f, μ)`](@ref) becomes equivalent to `v`. Instead of sampling `∫(L, v)` we can now sample `∫(L∘f, μ)`. The generates sample points `X_μ` can then be transformed to sample points `X_v = f.(X_μ)`.

Transformation functions `f = vartransform(v, μ)` support `f_inverse = InverseFunctions.inverse(f)` and `x_v, ladj = ChangesOfVariables.with_logabsdet_jacobian(f, x_μ)`.
