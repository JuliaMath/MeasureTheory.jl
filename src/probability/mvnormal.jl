
# Multivariate Normal distribution

import StatsFuns
using LinearAlgebra
export MvNormal
using Random
import Base

struct MvNormal{N, T, I, J} <: ParameterizedMeasure{N, T}
    cache::Vector{Float64}
    par::NamedTuple{N, T}
end

function getcache(d::MvNormal, n)
    view(getfield(d, :cache), 1:n)
end

function MvNormal(nt::NamedTuple{N,T}) where {N,T}
    I,J = mvnormaldims(nt)

    cache = Vector{Float64}(undef, max(I,J))
    MvNormal{N,T,I,J}(cache, nt)
end

function Base.size(d::MvNormal{N, T, I, J}) where {N,T,I,J}
    return (I,)
end

mvnormaldims(nt::NamedTuple{(:A, :b)}) = size(nt.A)

function MeasureTheory.basemeasure(μ::MvNormal{N, T, I,I}) where {N, T, I}
    return (1 / sqrt2π)^I * Lebesgue(ℝ)^I
end

sampletype(d::MvNormal{N, T, I, J}) where {N,T,I,J} = Vector{Float64}

MvNormal(; kwargs...) = begin
    MvNormal((; kwargs...))
end



function Random.rand!(rng::AbstractRNG, d::MvNormal{(:A, :b),T,I,J}, x::AbstractArray) where {T,I,J}
    z = getcache(d, J)
    rand!(rng, Normal()^J, z)
    copyto!(x, d.b)
    mul!(x, d.A, z, 1.0, 1.0)
    return x
end

function logdensity(d::MvNormal{(:A,:b)}, x)
    cache = getcache(d, size(d))
    z = d.A \ (x - d.b)
    return logdensity(MvNormal(), z) - logabsdet(d.A)
end

≪(::MvNormal, ::Lebesgue{ℝ}) = true
representative(::MvNormal{N, T, I,I}) where {N, T, I} = Lebesgue(ℝ)^I
