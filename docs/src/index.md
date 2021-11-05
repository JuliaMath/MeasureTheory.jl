```@meta
CurrentModule = MeasureTheory
```

# Home

`MeasureTheory.jl` is a package for building and reasoning about measures.

## Why?

A distribution (as provided by `Distributions.jl`) is also called a _probability measure_, and carries with it the constraint of adding (or integrating) to one. Statistical work usually requires this "at the end of the day", but enforcing it at each step of a computation can have considerable overhead. For instance, Bayesian modeling often requires working with unnormalized posterior densities or improper priors.

As a generalization of the concept of volume, measures also have applications outside of probability theory.

## Getting started

To install `MeasureTheory.jl`, open the Julia Pkg REPL (by typing `]` in the standard REPL) and run

```julia
pkg> add MeasureTheory
```

To get an idea of the possibilities offered by this package, go to the [documentation](https://cscherrer.github.io/MeasureTheory.jl/stable).

To know more about the underlying theory and its applications to probabilistic programming, check out our [JuliaCon 2021 submission](https://arxiv.org/abs/2110.00602).

## Support

[<img src=https://user-images.githubusercontent.com/1184449/140397787-9b7e3eb7-49cd-4c63-8f3c-e5cdc41e393d.png width="49%">](https://informativeprior.com/) [<img src=https://planting.space/sponsor/PlantingSpace-sponsor-3.png width=49%>](https://planting.space)