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
    m::T
end
kernel(m, ops...) = Kernel(ops, m)
kernel(m::T; ops...) where {T} = Kernel{T,typeof(ops.data)}(ops.data, m)
mapcall(t, x) = map(func -> func(x), t)
(k::Kernel)(x) = k.m(mapcall(k.ops, x)...)
(k::Kernel{<:Any,<:NamedTuple})(x) = k.m(;mapcall(k.ops, x)...)
