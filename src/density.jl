"""
    struct Density{M,B}
        μ::M
        base::B
    end

For measures μ and ν with μ≪ν, the density of μ with respect to ν (also called
the Radon-Nikodym derivative dμ/dν) is a function f defined on the support of ν
with the property that for any measurable a ⊂ supp(ν), μ(a) = ∫ₐ f dν.
    
Because this function is often difficult to express in closed form, there are
many different ways of computing it. We therefore provide a formal
representation to allow comptuational flexibilty.
"""
struct Density{M,B,L} <: Function
    μ::M
    base::B
    log::L
end

function 𝒹(μ::AbstractMeasure, ν::AbstractMeasure; log=true)
    return Density(μ, base; log=Val{log})
end

(f::Density{M,B,Val{true})(x) = logdensity(f.μ, x) - logdensity(f.base, x) 

(f::Density{M,B,Val{false})(x) = density(f.μ, x) / density(f, base(x))




"""
    struct DensityMeasure{F,B} <: AbstractMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct DensityMeasure{F,B,L} <: AbstractMeasure
    f    :: F
    base :: B
    log  :: L
end

∫(f, base::AbstractMeasure; log=true) = DensityMeasure(f, base; Val{log})

# TODO: `density` and `logdensity` functions for `DensityMeasure`
