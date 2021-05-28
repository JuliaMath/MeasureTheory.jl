
"""
    ≪(μ,ν)

# Absolute continuity

A measure μ is _absolutely continuous_ with respect to ν, written μ ≪ ν, if
ν(A)==0 implies μ(A)==0 for every ν-measurable set A.

Less formally, suppose we have a set A with ν(A)==0. If μ(A)≠0, then there can
be no way to "reweight" ν to get to μ. We can't make something from nothing.

This "reweighting" is really a density function. If μ≪ν, then there is some
function f that makes `μ == ∫(f,ν)` (see the help section for `∫`).

We can get this f directly via the Radon-Nikodym derivative, `f == 𝒹(μ,ν)` (see
the help section for `𝒹`).

Note that `≪` is not a partial order, because it is not antisymmetric. That is
to say, it's possible (in fact, common) to have two different measures `μ` and
`ν` with `μ ≪ ν` and `ν ≪ μ`. A simple example of this is 
```
μ = Normal()
ν = Lebesgue(ℝ)
```

When `≪` holds in both directions, the measures μ and ν are _equivalent_,
written `μ ≃ ν`. See the help section for `≃` for more information.
"""
function ≪ end


export ≃

"""
    ≃(μ,ν)

# Equivalence of Measure

Measures μ and ν on the same space X are equivalent, written `μ ≃ ν`, if `μ ≪ ν`
and `ν ≪ μ`. Note that this is often written `~` in the literature, but this is
overloaded in probabilistic programming, so we use this alternate notation. 

Also note that equivalence is very different from equality. For two equivalent
measures, the sets of non-zero measure will be identical, but what that measure
is in each case can be very different. 
"""
function ≃(μ,ν)
    return (μ≪ν && ν≪μ)
end

export representative

"""
    representative(μ::AbstractMeasure) -> AbstractMeasure

We need to be able to compute `μ ≪ ν` for each `μ` and `ν`. To do this directly
would require a huge number of methods (quadratic in the number of defined
measures). 

This function is a way around that. When defining a new measure `μ`, you should
also find some equivalent measure `ρ` that's "as primitive as possible". 

If possible, `ρ` should be a `PrimitiveMeasure`, or a `Product` of these. If
not, it should be a  transform (`Pushforward` or `Pullback`) of a
`PrimitiveMeasure` (or `Product` of these). 
"""
function representative(μ)
    function f(μ)
        # Check if we're done
        isprimitive(μ) && return μ
        ν = basemeasure(μ)
        return ν
    end

    fix(f, μ)
end

function ≪(μ, ν)
    μ == ν && return true
    representative(μ) ≪ representative(ν) && return true
    return false
end

@traitfn representative(μ::M) where {M; IsRepresentative{M}} = μ
