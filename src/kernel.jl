using KeywordCalls: instance_type


export Kernel

abstract type AbstractKernel <: AbstractMeasure end

struct Kernel{F,S} <: AbstractKernel
    f::F
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

kernel(μ, ops...) = Kernel(μ, ops)
kernel(μ, op) = Kernel(μ, op)

# kernel(Normal(μ=2))
function kernel(μ::P) where {P <: AbstractMeasure}
    (f, ops) = _kernelfactor(μ)
    Kernel{instance_type(f), typeof(ops)}(f,ops)
end

# kernel(Normal{(:μ,), Tuple{Int64}})
function kernel(::Type{P}) where {P <: AbstractMeasure}
    (f, ops) = _kernelfactor(P)
    Kernel{instance_type(f), typeof(ops)}(f,ops)
end



# kernel(::Type{P}, op::O) where {O, N, P<:ParameterizedMeasure{N}} = Kernel{constructor(P),O}(op)

# kernel(::Type{M}; ops...) where {M} = Kernel{Type{M},typeof(ops.data)}(M,ops.data)

# TODO: Would this benefit from https://github.com/tisztamo/FunctionWranglers.jl?
mapcall(t, x) = map(func -> func(x), t)

# (k::Kernel{Type{P},<:Tuple})(x) where {P<:ParameterizedMeasure} = k.f(mapcall(k.ops, x)...)

# (k::Kernel{M,<:NamedTuple})(x) where {M} = k.f(;mapcall(k.ops, x)...)

(k::Kernel)(x) = k.f(k.ops(x)...)

(k::Kernel{F,S})(x...) where {F, N, S<:NTuple{N,Symbol}} = k(x)


function (k::Kernel{F,S})(x::Tuple) where {F, N, S<:NTuple{N,Symbol}}
    k.f(NamedTuple{k.ops}(x))
end



# export kernelize

# function kernelize(μ::M) where {N, M <: ParameterizedMeasure{N}}
#     C = constructor(M)
#     (Kernel{C,}(NamedTuple{N}, ), values(getfield(μ, :par)))
# end

function _kernelfactor(::Type{P}) where {N, P <: ParameterizedMeasure{N}}
    (constructor(P), N)
end

function _kernelfactor(::P) where {N, P <: ParameterizedMeasure{N}}
    (constructor(P), N)
end

function _kernelfactor(μ::P) where {N, D<:ParameterizedMeasure{N}, P <: PowerMeasure{D}}
    C = constructor(D)
    (p -> C(p) ^ size(μ.data), N)
end
