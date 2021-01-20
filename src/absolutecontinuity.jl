
"""
    ≪(μ,ν)

# Absolute continuity

The following are equivalent:
1. μ ≪ ν
2. μ is absolutely continuous wrt ν
3. There exists a function f such that μ = ∫f dν
"""
function ≪ end

"""
    ≃(μ,ν)

# Equivalence of Measure

Measures μ and ν on the same space X are equivalent, written `μ ≃ ν`, if `μ ≪ ν` and `ν ≪ μ`.

Note that this is often written `~` in the literature, but this is overloaded in probabilistic programming, so we use alternate notation.
"""
function ≃(μ,ν)
    return (μ≪ν && ν≪μ)
end

export representative

function representative(μ)
    # Check if we're done
    isprimitive(μ) && return μ

    ν = basemeasure(μ)
    
    # Make sure not to leave the equivalence class
    (ν ≪ μ) || return μ

    # Fall back on a recusive call
    return representative(ν)
end

# TODO: ≪ needs more work
function ≪(μ, ν)
    μ == ν && return true
    representative(μ) ≪ representative(ν) && return true
end
