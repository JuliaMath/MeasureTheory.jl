using TransformVariables

export asparams

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

asparams(μ::ParameterizedMeasure, v::Val) = asparams(constructor(μ), v)
asparams(μ, s::Symbol) = asparams(μ, Val(s))

asparams(M::Type{A}) where {A<:AbstractMeasure} = asparams(M, NamedTuple())

function asparams(::Type{M}, constraints::NamedTuple{N}) where {N, M<: ParameterizedMeasure} 
    # @show M
    thekeys = params(M, constraints)
    t1 = NamedTuple{thekeys}(asparams(M, Val(k)) for k in thekeys)
    t2 = NamedTuple{N}(map(asConst, values(constraints)))
    C = constructorof(M)
    # @show C
    # @show constraints
    # @show transforms
    # Make sure we end up with a consistent ordering
    ordered_transforms = params(C(merge(t1, t2)))
    return TV.as(ordered_transforms)
end

# TODO: Make this work for Affine measures
# function asparams(::Type{M}, constraints::NamedTuple{N}) where {N, M<: Affine} 
#     ...
# end

export params

params(μ::ParameterizedMeasure) = getfield(μ, :par)

function params(::Type{M}, constraints::NamedTuple{N2}=NamedTuple()) where {N1, N2, M<: ParameterizedMeasure{N1}} 
thekeys = tuple((k for k in N1 if k ∉ N2)...)
end

params(μ) = ()

params(::Type{PM}) where {N, PM<:ParameterizedMeasure{N}} = N

params(::Type{PM}, ::NamedTuple{C}) where {N, C, PM<:ParameterizedMeasure{N}} = tuple(setdiff(N,C)...)

asparams(μ::ParameterizedMeasure, nt::NamedTuple=NamedTuple()) = asparams(constructor(μ), nt)

as(::Half) = asℝ₊

asparams(::Affine, ::Val{:μ}) = asℝ
asparams(::Affine, ::Val{:σ}) = asℝ₊
asparams(::Type{A}, ::Val{:μ}) where {A<:Affine} = asℝ
asparams(::Type{A}, ::Val{:σ}) where {A<:Affine} = asℝ₊

params(::Type{A}, nt::NamedTuple{C}) where {A<:Affine{N,M}} where {N,M,C} = tuple(setdiff(union(N, params(M)),C)...)
