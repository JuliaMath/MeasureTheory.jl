export Kernel

struct Kernel{T,S}
    ops::S
end

"""
    kernel(f, M)
    kernel((f1, f2, ...), M)

A kernel `κ = kernel(f, m)` returns a wrapper around
a function `f` giving the parameters for a measure of type `M`,
such that `κ(x) = M(f(x)...)`
respective `κ(x) = M(f1(x), f2(x), ...)`

If the argument is a named tuple `(;a=f1, b=f1)`, `κ(x)` is defined as
`M(;a=f(x),b=g(x))`.

# Reference

* https://en.wikipedia.org/wiki/Markov_kernel
"""
function kernel end

export kernel
kernel(op, ::Type{M}) where {M} = Kernel{M,typeof(op)}(op)
kernel(::Type{M}; ops...) where {M} = Kernel{M,typeof(ops.data)}(ops.data)
mapcall(t, x) = map(func -> func(x), t)
(k::Kernel{M,<:Tuple})(x) where {M} = M(mapcall(k.ops, x)...)
(k::Kernel{M,<:NamedTuple})(x) where {M} = M(;mapcall(k.ops, x)...)
(k::Kernel{M})(x) where {M} = M(k.ops(x)...)

export kernelize

function kernelize(μ::M) where {N, M <: ParameterizedMeasure{N}}
    C = M.name.wrapper
    (kernel(NamedTuple{N}, C), values(getfield(μ, :par)))
end
