"""
    kernel(f, m)
    kernel((f, g), m)

A kernel `κ = kernel(f, m)` returns a wrapper around
a function `f` giving the parameters for a measure `m`,
such that `κ(x) = m(f(x))`.
If `(a=f, b=g)` is a name tuple, `κ(x)` is defined as
`m(;a=f(x),b=g(x))`.

# Reference

* https://en.wikipedia.org/wiki/Markov_kernel
"""
struct Kernel{T,S}
    ops::S
end
kernel(::Type{M}, ops...) where {M} = Kernel{M,typeof(ops)}(ops)
kernel(::Type{M}; ops...) where {M} = Kernel{M,typeof(ops.data)}(ops.data)
mapcall(t, x) = map(func -> func(x), t)
(k::Kernel{M})(x) where {M} = M(mapcall(k.ops, x)...)
(k::Kernel{M,<:NamedTuple})(x) where {M} = M(;mapcall(k.ops, x)...)
