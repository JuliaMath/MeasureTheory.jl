
export Dirac

struct Dirac{X} <: AbstractMeasure
    x::X
end

sampletype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end
isprimitive(::Dirac) = true


Base.rand(_::Random.AbstractRNG, ::Type, μ::Dirac) = μ.x
