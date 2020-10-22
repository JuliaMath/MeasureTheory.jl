"""
    kernel(f, m)
    kernel((f, g), m)

A kernel `κ = kernel(f, m)` returns a wrapper around
a function `f` giving the parameters for a measure `m`,
such that `κ(x) = m(f(x))`.
If `(f, g)` is a tuple, `κ(x)` is defined as `merge`d return values.

# Reference

* https://en.wikipedia.org/wiki/Markov_kernel
"""
struct Kernel{T,S}
    ops::S
    m::T
end
kernel(m, ops...) = Kernel(ops, m)
kernel(m; ops...) = Kernel(ans.data, m)
mapcall(t, x) = map(func -> func(x), t)
(k::Kernel)(x) = k.m(mapcall(k.ops, x)...)
(k::Kernel{<:NamedTuple})(x) = k.m(;mapcall(k.ops, x)...)
