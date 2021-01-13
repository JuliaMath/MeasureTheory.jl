# Gumbel distribution

import StatsFuns
export Gumbel

@measure Gumbel(μ,σ) ≃ Lebesgue(ℝ)

function logdensity(d::Gumbel{()} , x::MaybeSym{T}) where {T <: Number}  
    return -exp(-x) - x
end
 
sampletype(::Gumbel{()}) = Real

import Base

function Base.rand(rng::AbstractRNG, d::Gumbel{()})
    u = rand()
    log(-log(u))
end

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true
representative(::Gumbel) = Lebesgue(ℝ)

@μσ_methods Gumbel()
