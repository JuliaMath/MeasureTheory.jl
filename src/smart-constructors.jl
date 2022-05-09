###############################################################################
# Affine

affine(f::AffineTransform, μ::AbstractMeasure) = Affine(f, μ)

affine(nt::NamedTuple, μ::AbstractMeasure) = affine(AffineTransform(nt), μ)

affine(f) = Base.Fix1(affine, f)

function affine(f::AffineTransform, parent::WeightedMeasure)
    WeightedMeasure(parent.logweight, affine(f, parent.base))
end

using MeasureBase: RealNumbers

function affine(f::AffineTransform{(:μ, :σ)}, parent::Lebesgue{RealNumbers})
    affine(AffineTransform((σ = f.σ,)), parent)
end

function affine(f::AffineTransform{(:μ, :ω)}, parent::Lebesgue{RealNumbers})
    affine(AffineTransform((ω = f.ω,)), parent)
end
