import Distributions

export Dists
const Dists = Distributions

@implement IsMeasure{D} where D <: Dists.Distribution

@implement HasDensity{D,X} where {VF, X, D <: Dists.Distribution{VF,Dists.Continuous}} begin
    baseMeasure(d) = Lebesgue(X)
    logdensity(d,x) = Dists.logpdf(d,x)
end

@implement HasRand{S} where {S <: Dists.Sampleable} begin
    rand(d) = Base.rand(d)
end
