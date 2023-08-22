###############################################################################
# AffinePushfwd

affine(f::AffineTransform, μ::AbstractMeasure) = AffinePushfwd(f, μ)

affine(nt::NamedTuple, μ::AbstractMeasure) = affine(AffineTransform(nt), μ)

affine(f) = Base.Fix1(affine, f)

function affine(f::AffineTransform, parent::WeightedMeasure)
    WeightedMeasure(parent.logweight, affine(f, parent.base))
end

using MeasureBase: RealNumbers

function affine(f::AffineTransform{(:μ, :σ)}, parent::Lebesgue{RealNumbers})
    affine(AffineTransform((σ = f.σ,)), parent)
end

function affine(f::AffineTransform{(:μ, :λ)}, parent::Lebesgue{RealNumbers})
    affine(AffineTransform((λ = f.λ,)), parent)
end

affine(f::AffineTransform, μ::LebesgueBase) = weightedmeasure(-logjac(f), μ)
