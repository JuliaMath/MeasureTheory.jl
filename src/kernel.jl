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
kernel(ops::Tuple, m) = Kernel(ops, m)
kernel(op, m) = Kernel((op,), m)

tuplemerge(x::Tuple) = x
tuplemerge(x...) = merge(x...)

(k::Kernel)(x) = k.m(tuplemerge((op(x) for op in k.ops)...)...)
