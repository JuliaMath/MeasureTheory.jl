
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

logdensity(μ::Dirac, x) = (x == μ.x) ? 0.0 : -Inf


Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x


export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))
