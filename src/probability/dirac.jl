
export Dirac

struct Dirac{X} <: AbstractMeasure
    supp :: X
end

Dirac(x::X) = Dirac{X}(x)

sampletype(μ::Dirac{X}) = X

function (μ::Dirac{X})(s) where {X}
    μ.supp ∈ s && return 1
    return 0
end

function logdensity(μ::Dirac{X}, x::X)
    μ.supp ∈ s && return 0.0
    return -Inf
end

Base.rand(μ::Dirac) = μ.supp
