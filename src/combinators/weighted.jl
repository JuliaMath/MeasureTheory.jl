export WeightedMeasure

"""
    struct WeightedMeasure{R,M} <: AbstractMeasure
        logweight :: R
        base :: M
    end
"""


abstract type AbstractWeightedMeasure <: AbstractMeasure end

logweight(μ::AbstractWeightedMeasure) = μ.logweight
basemeasure(μ::AbstractWeightedMeasure) = μ.base


# TV.as(μ::AbstractWeightedMeasure) = TV.as(μ.base)

function logdensity(sm::AbstractWeightedMeasure, x)
    logdensity(sm.base, x) + sm.logweight
end

###############################################################################

struct WeightedMeasure{R,M} <: AbstractWeightedMeasure
    logweight :: R
    base :: M
end

function Base.show(io::IO, ::MIME"text/plain", μ::WeightedMeasure)
    io = IOContext(io, :compact => true)
    print(io, exp(μ.logweight), " * ", μ.base)
end

function Base.show_unquoted(io::IO, μ::WeightedMeasure, indent::Int, prec::Int)
    io = IOContext(io, :compact => true)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end

function Base.:*(k::T, m::AbstractMeasure) where {T <: Number}
    logk = log(k)
    return WeightedMeasure{typeof(logk), typeof(m)}(logk,m)
end

Base.:*(m::AbstractMeasure, k::Real) = k * m

≪(::M, ::WeightedMeasure{R,M}) where {R,M} = true
≪(::WeightedMeasure{R,M}, ::M) where {R,M} = true





sampletype(μ:: WeightedMeasure) = sampletype(μ.base)

###############################################################################
struct ParamWeightedMeasure{F,N,T,R,B} <: AbstractWeightedMeasure 
    f::F
    par::NamedTuple{N,T}
    logweight::R
    base::B

    function ParamWeightedMeasure(f::F, par::NamedTuple{N,T}, base::B) where {F,N,T,B}
        ℓ = f(par)
        new{F,N,T,typeof(ℓ),B}(f,par,ℓ,base)
    end
end

function Base.show(io::IO, ::MIME"text/plain", d::ParamWeightedMeasure)
    io = IOContext(io, :compact => true)
    print(io, "ParamWeighted(",d.f, ", ", d.par,", ", d.base, ")")
end
