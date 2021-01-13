
# Cauchy distribution

import StatsFuns
export Cauchy

import Base: eltype


@measure Cauchy(μ,σ) ≃ (1/π) * Lebesgue(ℝ)


# Standard Cauchy

function logdensity(d::Cauchy{()} , x::X) where {X}  
    return -log(1 + x^2)
end

# Cauchy() = Cauchy{EmptyNamedTuple,Real}(NamedTuple())
 
sampletype(::Cauchy) = Real

Base.rand(μ::Cauchy{()}) = rand(Dists.Cauchy())

≪(::Cauchy, ::Lebesgue{X}) where X <: Real = true
representative(::Cauchy) = Lebesgue(ℝ)

@μσ_methods Cauchy()
