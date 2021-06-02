export PrimitiveMeasure

abstract type PrimitiveMeasure <: AbstractMeasure end

basemeasure(μ::PrimitiveMeasure) = μ

logdensity(μ::PrimitiveMeasure, x) = 0.0

logdensity(μ::M, ν::M, x) where {M<:PrimitiveMeasure} = 0.0
