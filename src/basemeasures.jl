const BASE_MEASURES = SimplePoset()

export baseMeasure

function baseMeasure

"""
    ≪(μ,ν)

# Absolute continuity

The following are equivalent:
1. μ ≪ ν
2. μ is absolutely continuous wrt ν
3. There exists a function f such that μ = ∫f dν
"""
function ≪(μ,ν)
    μ == ν && return true
    
    a = baseMeasure(μ)
    b = baseMeasure(ν)
    a == b && return true

    return b ∈ BASE_MEASURES.D.N[a]
end

"""
    ≅(μ,ν)

# Equivalence of Measure

Measures μ and ν on the same space X are equivalent, written `μ ≅ ν`, if `μ ≪ ν` and `ν ≪ μ`.

Note that this is often written `~` in the literature, but this is overloaded in probabilistic programming, so we use alternate notation.
"""
function ≅(μ,ν)
    return (μ≪ν && ν≪μ)
end
