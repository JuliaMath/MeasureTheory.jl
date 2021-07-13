using TransformVariables

export ParameterizedMeasure
abstract type ParameterizedMeasure{N} <: AbstractMeasure end

function Base.getproperty(μ::ParameterizedMeasure{N}, prop::Symbol) where {N}
    return getproperty(getfield(μ, :par), prop)
end

function Base.propertynames(μ::ParameterizedMeasure{N}) where {N}
    return N
end

function Base.show(io::IO, ::MIME"text/plain", μ::ParameterizedMeasure{()}) 
    io = IOContext(io, :compact => true)
    print(io, nameof(typeof(μ)), "()")
end

function Base.show(io::IO, ::MIME"text/plain", μ::ParameterizedMeasure{N}) where {N}
    io = IOContext(io, :compact => true)
    print(io, nameof(typeof(μ)))
    print(io, getfield(μ,:par))
end

export asparams

# Allow things like
#
# julia> Normal{(:μ,)}(2)
# Normal(μ = 2,)
#
function (::Type{P})(args...) where {N, P <: ParameterizedMeasure{N}}
    C = constructorof(P)
    return C(NamedTuple{N}(args...))
end

(::Type{P})(;kwargs...) where {P <: ParameterizedMeasure} = P(NamedTuple(kwargs))

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

asparams(μ::ParameterizedMeasure, v::Val) = asparams(typeof(μ), v)
asparams(μ, s::Symbol) = asparams(μ, Val(s))

asparams(M::Type{PM}) where {PM<:ParameterizedMeasure} = asparams(M, NamedTuple())

function asparams(::Type{M}, constraints::NamedTuple{N2}) where {N1, N2, M<: ParameterizedMeasure{N1}} 
    # @show M
    thekeys = params(M, constraints)
    transforms = NamedTuple{thekeys}(asparams(M, Val(k)) for k in thekeys)
    C = constructorof(M)
    # @show C
    # @show constraints
    # @show transforms
    # Make sure we end up with a consistent ordering
    ordered_transforms = NamedTuple{thekeys}(params(C(merge(constraints, transforms))))
    return as(ordered_transforms)
end


asparams(μ::ParameterizedMeasure, nt::NamedTuple=NamedTuple()) = asparams(typeof(μ), nt)

export params

params(μ::ParameterizedMeasure) = getfield(μ, :par)

function params(::Type{M}, constraints::NamedTuple{N2}=NamedTuple()) where {N1, N2, M<: ParameterizedMeasure{N1}} 
thekeys = tuple((k for k in N1 if k ∉ N2)...)
end

params(μ) = ()


function ConstructionBase.setproperties(d::P, nt::NamedTuple) where {P<:ParameterizedMeasure}
    return constructorof(P)(merge(params(d), nt)) 
end

function Accessors.set(d::ParameterizedMeasure, ::typeof(params), nt::NamedTuple)
    setproperties(d, nt)
end

function Accessors.set(d::ParameterizedMeasure{N}, ::typeof(params), p) where {N}
    setproperties(d, NamedTuple{N}(p...))
end
