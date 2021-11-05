```@meta
CurrentModule = MeasureTheory
```

# MeasureTheory

## Introduction

There are lots of packages for working with probability distributions. But very
often, we need to work with "distributions" that really aren't. 

For example, the [correspondence between regularization and Bayesian prior
distributions](https://en.wikipedia.org/wiki/Bayesian_interpretation_of_kernel_regularization)
leads naturally to the idea of extending probabilistic programming systems to
cover both. But it's easy to come up with a loss function for which the integral
of the
corresponding "prior" is infinite! The result is not really a distirbution. It
is, however, still a measure.

Even restricted to Bayesian methods, users might sometimes want to use an
[improper
prior](https://en.wikipedia.org/wiki/Prior_probability#Improper_priors). By
definition, these cannot be integrated over their domain. But an improper prior
is still a measure.

In [Markov chain Monte Carlo
(MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo), we often work
with distributions for which we can only caluculate 
the log-density up to an additive constant. Considering this instead as a
measure can be helpful. Even better, consdering intermediate computations along
the way as computations on measures saves us from computing normalization terms
where the end result will discard this anyway.

To be clear, that's not to say that we always discard normalizations. Rather,
they're considered as belonging to the measure itself, rather than being
included in each sub-computation. If measures you work with happen to also be
probability distributions, you'll always be able to recover those results.


## Index

```@index
```