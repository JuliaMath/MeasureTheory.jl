export WeightedMeasure

"""
    struct WeightedMeasure{R,M} <: AbstractMeasure
        logweight :: R
        base :: M
    end


"""
struct WeightedMeasure{R,M} <: AbstractMeasure
    logweight :: R
    base :: M
end

logweight(μ::WeightedMeasure) = μ.logweight

function Base.show(io::IO, μ::WeightedMeasure)
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

function logdensity(sm::WeightedMeasure{R,M}, x::X) where {X, R, M <: AbstractMeasure}
    logdensity(sm.base, x) + sm.logweight
end

function Base.:*(k::T, m::AbstractMeasure) where {T <: Number}
    logk = log(k)
    return WeightedMeasure{typeof(logk), typeof(m)}(logk,m)
end

Base.:*(m::AbstractMeasure, k::Real) = k * m

≪(::M, ::WeightedMeasure{R,M}) where {R,M} = true
≪(::WeightedMeasure{R,M}, ::M) where {R,M} = true

basemeasure(μ::WeightedMeasure{R,M}) where {R,M} = μ.base

sampletype(μ:: WeightedMeasure) = sampletype(μ.base)
