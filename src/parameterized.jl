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

# TODO: Make this work for AffinePushfwd measures
# function asparams(::Type{M}, constraints::NamedTuple{N}) where {N, M<: AffinePushfwd} 
#     ...
# end

function asparams(μ::M, nt::NamedTuple = NamedTuple()) where {M<:ParameterizedMeasure}
    asparams(constructorof(M), nt)
end

as(::Half) = asℝ₊


asparams(::AffinePushfwd, ::StaticSymbol{:μ}) = asℝ
asparams(::AffinePushfwd, ::StaticSymbol{:σ}) = asℝ₊
asparams(::Type{A}, ::StaticSymbol{:μ}) where {A<:AffinePushfwd} = asℝ
asparams(::Type{A}, ::StaticSymbol{:σ}) where {A<:AffinePushfwd} = asℝ₊

function asparams(d::AffinePushfwd{N,M,T}, ::StaticSymbol{:μ}) where {N,M,T<:AbstractArray}
    as(Array, asℝ, size(d.μ))
end

function asparams(d::AffinePushfwd{N,M,T}, ::StaticSymbol{:σ}) where {N,M,T<:AbstractArray}
    as(Array, asℝ, size(d.σ))
end


asparams(::Type{<:Bernoulli}, ::StaticSymbol{:p}) = as𝕀
asparams(::Type{<:Bernoulli}, ::StaticSymbol{:logitp}) = asℝ

asparams(::Type{<:Beta}, ::StaticSymbol{:α}) = asℝ₊
asparams(::Type{<:Beta}, ::StaticSymbol{:β}) = asℝ₊


asparams(::Type{<:BetaBinomial}, ::StaticSymbol{:α}) = asℝ₊
asparams(::Type{<:BetaBinomial}, ::StaticSymbol{:β}) = asℝ₊


asparams(::Type{<:Binomial}, ::StaticSymbol{:p}) = as𝕀
asparams(::Type{<:Binomial}, ::StaticSymbol{:logitp}) = asℝ
asparams(::Type{<:Binomial}, ::StaticSymbol{:probitp}) = asℝ

asparams(::Type{<:Exponential}, ::StaticSymbol{:β}) = asℝ₊
asparams(::Type{<:Exponential}, ::StaticSymbol{:logβ}) = asℝ
asparams(::Type{<:Exponential}, ::StaticSymbol{:λ}) = asℝ₊
asparams(::Type{<:Exponential}, ::StaticSymbol{:logλ}) = asℝ


asparams(::Type{<:LKJCholesky}, ::StaticSymbol{:η}) = asℝ₊
asparams(::Type{<:LKJCholesky}, ::StaticSymbol{:logη}) = asℝ

asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:p}) = as𝕀
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:logitp}) = asℝ
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:r}) = asℝ₊
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:λ}) = asℝ₊
asparams(::Type{<:NegativeBinomial}, ::StaticSymbol{:logλ}) = asℝ

asparams(::Type{<:Normal}, ::StaticSymbol{:σ²}) = asℝ₊
asparams(::Type{<:Normal}, ::StaticSymbol{:τ}) = asℝ₊
asparams(::Type{<:Normal}, ::StaticSymbol{:logτ}) = asℝ

asparams(::Type{<:Poisson}, ::StaticSymbol{:λ}) = asℝ₊
asparams(::Type{<:Poisson}, ::StaticSymbol{:logλ}) = asℝ

asparams(::Type{<:SnedecorF}, ::StaticSymbol{:ν1}) = asℝ₊
asparams(::Type{<:SnedecorF}, ::StaticSymbol{:ν2}) = asℝ₊

asparams(::Type{<:StudentT}, ::StaticSymbol{:ν}) = asℝ₊