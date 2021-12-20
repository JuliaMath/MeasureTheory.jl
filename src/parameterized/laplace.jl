
# Laplace distribution

export Laplace

@parameterized Laplace() 

for N in AFFINEPARS
    @eval begin
        proxy(d::Laplace{$N}) = affine(params(d), Laplace())
    end
end

@inline function logdensity_def(d::Laplace{()}, x)
    return -abs(x)
end

logdensity_def(d::Laplace, x) = logdensity_def(proxy(d), x)


basemeasure(::Laplace{()}) = WeightedMeasure(static(-logtwo), Lebesgue(ℝ))
basemeasure(d::Laplace) = basemeasure(proxy(d))

function tbasemeasure_type(::Type{<:Laplace{()}}) 
    WeightedMeasure{Float64, Lebesgue{MeasureBase.RealNumbers}}
end


# @affinepars Laplace


Base.rand(rng::AbstractRNG, ::Type{T}, μ::Laplace{()}) where {T} =
    rand(rng, Dists.Laplace())
Base.rand(rng::AbstractRNG, ::Type{T}, μ::Laplace) where {T} = Base.rand(rng, T, proxy(μ))

≪(::Laplace, ::Lebesgue{X}) where {X<:Real} = true

TV.as(::Laplace) = asℝ
