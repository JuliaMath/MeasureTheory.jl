
export Dirac

struct Dirac{X} <: Measure{X}
    supp :: X
end

Dirac(x::X) = Dirac{X}(x)

function (μ::Dirac{X})(s) where {X}
    μ.supp ∈ s && return 1
    return 0
end
