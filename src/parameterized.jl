export asparams

"""
`asparams` build on `TransformVariables.as` to construct bijections to the
parameter space of a given parameterized measure. Because this is only possible
for continuous parameter spaces, we allow constraints to assign values to any
subset of the parameters.

--------

    asparams(::Type{<:ParameterizedMeasure}, ::StaticSymbol)

Return a transformation for a particular parameter of a given parameterized
measure. For example,

```
julia> asparams(Normal, static(:σ))
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

function asparams(μ::M, v::StaticSymbol) where {M<:ParameterizedMeasure}
    asparams(constructorof(M), v)
end

asparams(μ, s::Symbol) = asparams(μ, static(s))

asparams(M::Type{A}) where {A<:AbstractMeasure} = asparams(M, NamedTuple())

function asparams(::Type{M}, constraints::NamedTuple{N}) where {N,M<:ParameterizedMeasure}
    # @show M
    thekeys = paramnames(M, constraints)
    t1 = NamedTuple{thekeys}(asparams(M, StaticSymbol(k)) for k in thekeys)
    t2 = NamedTuple{N}(map(asConst, values(constraints)))
    C = constructorof(M)
    # @show C
    # @show constraints
    # @show transforms
    # Make sure we end up with a consistent ordering
    ordered_transforms = params(C(merge(t1, t2)))
    return as(ordered_transforms)
end

# TODO: Make this work for Affine measures
# function asparams(::Type{M}, constraints::NamedTuple{N}) where {N, M<: Affine} 
#     ...
# end

function asparams(μ::M, nt::NamedTuple = NamedTuple()) where {M<:ParameterizedMeasure}
    asparams(constructorof(M), nt)
end

as(::Half) = asℝ₊
