export asConst
struct AsConst{T} <: TV.VectorTransform
    value::T
end

asConst(x) = AsConst(x)

as(d::Dirac) = AsConst(d.x)

TV.dimension(t::AsConst) = 0

function TV.transform_with(flag::TV.NoLogJac, t::AsConst, x, index)
    return (t.value, TV.NoLogJac(), index)
end

function TV.transform_with(flag::TV.LogJac, t::AsConst, x, index)
    return (t.value, 0.0, index)
end

TV.inverse_eltype(t::AsConst, x::AbstractVector) = Any

function TV.inverse_at!(x::AbstractVector, index, t::AsConst, y::AbstractVector)
    return index
end
