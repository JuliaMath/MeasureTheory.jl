struct OrthoTransform{I,N,T} <: TV.AbstractTransform
    dim::I
    par::NamedTuple{N,T}
end

TV.dimension(t::OrthoTransform) = t.dim

TV.inverse_eltype(t::OrthoTransform, y::AbstractVector) = extended_eltype(y)

function TV.inverse_at!(
    x::AbstractVector,
    index,
    t::OrthoTransform{(:σ,)},
    y::AbstractVector,
) end

function TV.transform_with(
    flag::LogJacFlag,
    t::OrthoTransform{(:σ,)},
    x::AbstractVector,
    index,
)
    ylen, xlen = size(d.σ)
    T = extended_eltype(x)
    ℓ = logjac_zero(flag, T)
    y = d.σ * view(x, index, index + xlen - 1)
    index += xlen
    y, ℓ, index
end
