using LinearAlgebra: LinearAlgebra

"""
    logmass(μ::AbstractMeasure)

Compute the logarithm of the [`mass`](@ref) of the measure ``μ``.
"""
function logmass end

"""
    mass(μ::AbstractMeasure)

Compute the mass of the measure ``μ``.

That is, given a set ``S`` with elements ``x ∈ S``, compute
``μ(S) = ∫_S dμ(x)``.

!!! note
    This function should not be directly modified. Instead, developers should
    overload [`logmass`](@ref).
"""
mass(μ) = Exp(logmass(μ))

"""
    normalize(μ::AbstractMeasure)

Convert the measure ``μ`` to a probability measure, that is, a measure whose
total mass is 1.
"""
LinearAlgebra.normalize(μ::AbstractMeasure) = μ * inv(mass(μ))

"""
    is_probability_measure(μ::AbstractMeasure; kwargs...) -> Bool

Return `true` if ``μ`` is a probability measure, that is, a measure whose
total [`mass`](@ref) is 1.

`kwargs` are forwarded to `isapprox`.
"""
is_probability_measure(μ; kwargs...) = isapprox(mass(μ), true; kwargs...)
