export SpikeMixture

struct SpikeMixture{T,S} <: AbstractMeasure
    m::T
    w::S
    w̄::S
end
SpikeMixture(m, w) = SpikeMixture(m, w, 1-w)

function Base.show(io::IO, μ::SpikeMixture)
    io = IOContext(io, :compact => true)
    print(io, "(", μ.w, "*", string(μ.m), " + ", μ.w̄, "Dirac(0))")
end

basemeasure(μ::SpikeMixture) = SpikeMixture(basemeasure(μ.m), 1.0, (1/μ.w - 1)/density(μ.m, 0))
function MeasureTheory.logdensity(μ::SpikeMixture, x)
    return log(μ.w) + logdensity(μ.m, x)
end
function sampletype(μ::SpikeMixture)
    sampletype(μ.m)
end
function Base.rand(rng::AbstractRNG, T::Type, μ::SpikeMixture)
     return T((rand(rng) < μ.w)*rand(rng, μ.m))
end