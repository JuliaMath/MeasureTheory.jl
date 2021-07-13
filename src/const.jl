using TransformVariables
const TV = TransformVariables

struct Const{T} <: Function
    value::T
end

export asConst
struct AsConst{T} <: TV.VectorTransform
    value::T
end

asConst(x) = AsConst(x)

TV.as(c::Const) = AsConst(c.value)

function Base.show(io::IO, ::MIME"text/plain", c::Const)
    io = IOContext(io, :compact => true)
    print(io, "_ -> ")
    print(io, c.t)
end

(c::Const)(x) = c.value

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
