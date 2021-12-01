###############################################################################
# Affine

affine(f::AffineTransform, μ::AbstractMeasure) = Affine(f, μ)

affine(nt::NamedTuple, μ::AbstractMeasure) = affine(AffineTransform(nt), μ)

affine(f) = μ -> affine(f, μ)

function affine(f::AffineTransform, parent::WeightedMeasure)
    WeightedMeasure(parent.logweight, affine(f, parent.base))
end

function affine(f::AffineTransform, parent::FactoredBase)
    constℓ = parent.constℓ
    varℓ = parent.varℓ
    # Avoid transforming `inbounds`, which is expensive
    base = affine(f, restrict(parent.inbounds, parent.base))
    FactoredBase(Returns(true), constℓ, varℓ, base)
end
