
# Multivariate Normal distribution

using LinearAlgebra
export MvNormal
using Random
import Base


@parameterized MvNormal(μ, Σ)

@kwstruct MvNormal(Σ)

@kwstruct MvNormal(μ, Σ)

@kwstruct MvNormal(σ)

@kwstruct MvNormal(μ, σ)

@kwstruct MvNormal(μ, Ω)

@kwstruct MvNormal(Ω)

@kwstruct MvNormal(ω)
