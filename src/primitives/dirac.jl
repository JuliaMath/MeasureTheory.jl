
export Dirac

struct Dirac{X} <: PrimitiveMeasure
    x::X
end

sampletype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

logdensity(μ::Dirac, x) = (x == μ.x) ? 0.0 : -Inf


Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x


export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

function logdensity(μ::Dirac{M}, ν::Dirac{M}, x) where {M} 
    if μ.x == ν.x
        x == μ.x && return 0.0
        return -Inf
    elseif μ.x == x
        return Inf
    else
        return -Inf
    end
end

testvalue(d::Dirac) = d.x
