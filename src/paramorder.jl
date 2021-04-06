
# This list defines the order for parameters in `ParameterizedMeasure`s. If
# ordering for a given measure seems strange, it can be corrected by adjusting
# this list accordingly.
PARAMS = [
    :μ
    :σ
    :logσ
    :τ
    :Σ
    :σ²
    :n
    :r
    :p
    :logitp
    :probitp
    :λ
]

function getorder(mylist)
    result = Dict{Symbol, Float64}()
    for (j,par) in enumerate(mylist)
        result[par] = j
    end
    result
end

# NOTE: External packages can add new parameter names by assigning new entries
# to PARAM_ORDER. Values in this dictionary are Float64, to be sure there are
# plenty of "in-between" values.
const PARAM_ORDER = getorder(PARAMS)

using Accessors
using NamedTupleTools: namedtuple

_paramsort(nt::NamedTuple{(), Tuple{}}) = nt



function _paramsort(nt::NamedTuple{K,V}) where {K,V}
    # Assign each symbol a default rank
    namerank(p) = Base.Math.evalpoly(256, reverse(Float64.(collect(String(p)))))
    π = sortperm(collect(K), by = p -> get(PARAM_ORDER, p, namerank(p)))
    k = K[π]
    v = @inbounds (_paramsort.(values(nt)))[π]
    return namedtuple(k)(v)
end

_paramsort(t::Tuple) = _paramsort.(t)
_paramsort(x) = x

struct Lenses{T}
    ls::T
end

call(f,args...; kwargs...) = f(args...; kwargs...)

(ℓ::Lenses)(nt::NamedTuple) = Base.Fix2(call, nt).(ℓ.ls)

@generated function paramsort(nt::NamedTuple{K,V}) where {K,V}
    s = _paramsort(schema(nt))
    ℓ = Lenses(lenses(s))
    return :(leaf_setter($s)($ℓ(nt)...))
end

export paramsort

paramsort(t::T) where {T<:Tuple} = paramsort.(t)
paramsort(x) = x

function paramsort(lm::NestedTuples.LazyMerge)
    x = getfield(lm, :x)
    y = getfield(lm, :y)
    return lazymerge(paramsort(x), paramsort(y))
end
