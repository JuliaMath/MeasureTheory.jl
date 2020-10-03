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
struct Density{M,B} <: Function
    μ::M
    base::B
end

"""
    struct DensityMeasure{F,B} <: AbstractMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct DensityMeasure{X,F,B} <: AbstractMeasure
    density :: F
    base    :: B
end

# function density(μ::M, ν::M) where {M}
#     if  μ==ν
#         return () -> 1.0
#     end
# end

density(μ, base=baseMeasure(μ)) = Density(μ, base)
logdensity(μ, base=baseMeasure(μ)) = LogDensity(μ, base)

density(μ::Dists.Distribution{Dists.Univariate,Dists.Continuous}, x::Real) = pdf(μ,x)
logdensity(μ::Dists.Distribution{Dists.Univariate,Dists.Continuous}, x::Real) = logpdf(μ,x)

density(μ::AbstractMeasure, x::X) where {X} = density(μ, baseMeasure(μ))(x) 

logdensity(μ::AbstractMeasure, y::Y) where {X, Y <: X} = logdensity(μ, baseMeasure(μ))(x)
