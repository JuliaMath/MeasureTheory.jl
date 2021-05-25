using TransformVariables

export ParameterizedMeasure
abstract type ParameterizedMeasure{N} <: AbstractMeasure end

function Base.getproperty(μ::ParameterizedMeasure{N}, prop::Symbol) where {N}
    return getproperty(getfield(μ, :par), prop)
end

function Base.propertynames(μ::ParameterizedMeasure{N}) where {N}
    return N
end

function Base.show(io::IO, μ::ParameterizedMeasure{()}) 
    print(io, nameof(typeof(μ)), "()")
end

function Base.show(io::IO, μ::ParameterizedMeasure{N}) where {N}
    io = IOContext(io, :compact => true)
    print(io, nameof(typeof(μ)))
    print(io, getfield(μ,:par))
end

# e.g. Normal(;μ=μ, σ=σ) = Normal((μ=μ, σ=σ))
(M::Type{<: ParameterizedMeasure})(; kwargs...) = M(NamedTuple(kwargs))

(M::Type{<: ParameterizedMeasure})(::Tuple{}) = M(NamedTuple())

export asparams

# Allow things like
#
# julia> Normal{(:μ,)}(2)
# Normal(μ = 2,)
#
function (M::Type{P})(args...) where {N, P <: ParameterizedMeasure{N}}
    constructor = M.body.name.wrapper
    return constructor(NamedTuple{N}(args...))
end

"""
`asparams` build on `TransformVariables.as` to construct bijections to the
parameter space of a given parameterized measure. Because this is only possible
for continuous parameter spaces, we allow constraints to assign values to any
subset of the parameters.

--------

    asparams(::Type{<:ParameterizedMeasure}, ::Val{::Symbol})

Return a transformation for a particular parameter of a given parameterized
measure. For example,

```
julia> asparams(Normal, Val(:σ))
asℝ₊
```

-----

    asparams(::Type{<: ParameterizedMeasure{N}}, constraints::NamedTuple) where {N}

Return a transformation for a given parameterized measure subject to the named tuple
`constraints`. For example,

```
julia> asparams(Binomial{(:p,)}, (n=10,))
TransformVariables.TransformTuple{NamedTuple{(:p,), Tuple{TransformVariables.ScaledShiftedLogistic{Float64}}}}((p = as𝕀,), 1)
```

------------

    aspararams(::ParameterizedMeasure)

Return a transformation with no constraints. For example,

```
julia> asparams(Normal{(:μ,:σ)})
TransformVariables.TransformTuple{NamedTuple{(:μ, :σ), Tuple{TransformVariables.Identity, TransformVariables.ShiftedExp{true, Float64}}}}((μ = asℝ, σ = asℝ₊), 2)
```
"""
function asparams end

asparams(M::Type{PM}) where {PM<:ParameterizedMeasure} = asparams(M, NamedTuple())

function asparams(M::Type{<: ParameterizedMeasure{N}}, constraints::NamedTuple) where {N} 
    thekeys = setdiff(N, keys(constraints))

    result = NamedTuple()
    for k in thekeys
        t = asparams(M, Val(k))
        result = merge(result, NamedTuple{(k,)}((t,)))
    end
    return as(paramsort(result))
end

export params

params(μ::ParameterizedMeasure) = getfield(μ, :par)
