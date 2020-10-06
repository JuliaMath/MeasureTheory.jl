"""
    struct ScaledMeasure{R,M} <: AbstractMeasure
        logscale :: R
        base :: M
    end
    

"""
struct ScaledMeasure{R,M} <: AbstractMeasure
    logscale :: R
    base :: M
end

Base.show(io::IO, μ::ScaledMeasure) = print(io, exp(μ.logscale), " * ", μ.base)

function Base.show_unquoted(io::IO, μ::ProductMeasure, indent::Int, prec::Int)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end

function logdensity(sm::ScaledMeasure{R,M}, x::X) where {X, R, M <: AbstractMeasure}
    logdensity(sm.base, x) + sm.logscale
end

function Base.:*(k::Number, m::AbstractMeasure) 
    if iszero(k)
        return TrivialMeasure
    else
        return ScaledMeasure{typeof(k), typeof(m)}(log(k),m)
    end
end

Base.:*(m::AbstractMeasure, k::Real) = k * m

≪(::M, ::ScaledMeasure{R,M}) where {R,M} = true
≪(::ScaledMeasure{R,M}, ::M) where {R,M} = true

basemeasure(μ::ScaledMeasure{R,M}) where {R,M} = μ.base

sampletype(μ:: ScaledMeasure) = sampletype(μ.base)
