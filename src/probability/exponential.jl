
# Exponential distribution

import StatsFuns
export Exponential

@measure Exponential(λ) ≪ Lebesgue(ℝ₊)

function logdensity(d::Exponential{()} , x)
    return -x
end

sampletype(::Exponential) = Real

Base.rand(μ::Exponential{()}) = -log(rand())

≪(::Exponential, ::Lebesgue{ℝ₊}) = true

representative(::Exponential) = Lebesgue(ℝ₊)
