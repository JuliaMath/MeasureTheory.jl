# MeasureTheory

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cscherrer.github.io/MeasureTheory.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cscherrer.github.io/MeasureTheory.jl/dev)
[![Build Status](https://github.com/cscherrer/MeasureTheory.jl/workflows/CI/badge.svg)](https://github.com/cscherrer/MeasureTheory.jl/actions)
[![Coverage](https://codecov.io/gh/cscherrer/MeasureTheory.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cscherrer/MeasureTheory.jl)

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

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

## Stargazers over time

[![Stargazers over time](https://starchart.cc/cscherrer/MeasureTheory.jl.svg)](https://starchart.cc/cscherrer/MeasureTheory.jl)
