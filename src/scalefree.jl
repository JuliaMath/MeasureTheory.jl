export ScaleFree    
using MeasureBase: LebesgueBase
struct ScaleFree <: AbstractMeasure end

MeasureBase.insupport(::ScaleFree, x) = x > zero(x)

@inline MeasureBase.logdensity_def(::ScaleFree, x) = inv(x)
@inline MeasureBase.basemeasure(::ScaleFree) = LebesgueBase()

MeasureBase.testvalue(::Type{T}, d::ScaleFree) where {T} = one(T)

function Base.rand(rng::AbstractRNG, ::Type{T}, d::ScaleFree) where {T}
    one(T)
end
