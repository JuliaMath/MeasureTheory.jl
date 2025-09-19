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
TransformVariables.TransformTuple{NamedTuple{(:p,), Tuple{TransformVariables.CompositeScalarTransform{Tuple{TVShift{Float64}, TVScale{Float64}, TVLogistic}}}}}((p = as𝕀,), 1)
```

------------

    aspararams(::ParameterizedMeasure)

Return a transformation with no constraints. For example,

```
julia> asparams(Normal{(:μ,:σ)})
TransformVariables.TransformTuple{NamedTuple{(:μ, :σ), Tuple{TransformVariables.Identity, TransformVariables.CompositeScalarTransform{Tuple{TVShift{Float64}, TVExp}}}}}((μ = asℝ, σ = asℝ₊), 2)
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

asparams(::Type{<:Normal}, ::StaticSymbol{:μ}) = asℝ
asparams(::Type{<:Normal}, ::StaticSymbol{:σ}) = asℝ₊
asparams(::Type{<:Normal}, ::StaticSymbol{:σ²}) = asℝ₊
asparams(::Type{<:Normal}, ::StaticSymbol{:τ}) = asℝ₊
asparams(::Type{<:Normal}, ::StaticSymbol{:logτ}) = asℝ

asparams(::Type{<:Poisson}, ::StaticSymbol{:λ}) = asℝ₊
asparams(::Type{<:Poisson}, ::StaticSymbol{:logλ}) = asℝ

asparams(::Type{<:SnedecorF}, ::StaticSymbol{:ν1}) = asℝ₊
asparams(::Type{<:SnedecorF}, ::StaticSymbol{:ν2}) = asℝ₊

asparams(::Type{<:StudentT}, ::StaticSymbol{:ν}) = asℝ₊

function as(d::PowerMeasure)
    as(Array, as(d.parent), length.(d.axes)...)
end

function as(d::ProductMeasure{<:AbstractArray{<:Dirac}})
    return asConst(testvalue.(marginals(d)))
end

function as(d::ProductMeasure{A}) where {A<:AbstractArray}
    mar = marginals(d)
    ts = map(as, mar)
    if allequal(ts)
        return as(Array, first(ts), size(ts))
    else
        error("Not yet implemented")
    end
end

function as(d::ProductMeasure{A}) where {A<:MappedArrays.ReadonlyMappedArray}
    d1 = marginals(d).f(first(marginals(d).data))
    as(Array, as(d1), size(marginals(d))...)
end

function as(d::ProductMeasure{T}) where {T<:Tuple}
    as(map(as, d.marginals))
end

function as(d::ProductMeasure{<:Base.Generator})
    d1 = marginals(d).f(first(marginals(d).iter))
    as(Array, as(d1), size(marginals(d))...)
end

as(::Beta) = as𝕀
as(::Cauchy) = asℝ
as(d::Dirichlet{(:α,)}) = TV.UnitSimplex(length(d.α))
as(::Exponential) = asℝ₊
as(::Gamma) = asℝ₊
as(::Gumbel) = asℝ
as(::InverseGaussian) = asℝ₊
as(::Laplace) = asℝ
as(d::MvNormal{(:μ,)}) = as(Array, length(d.μ))

as(d::MvNormal{(:Σ,),Tuple{C}}) where {C<:Cholesky} = as(Array, size(d.Σ, 1))
as(d::MvNormal{(:Λ,),Tuple{C}}) where {C<:Cholesky} = as(Array, size(d.Λ, 1))
as(d::MvNormal{(:μ, :Σ),<:Tuple{T,C}}) where {T,C<:Cholesky} = as(Array, size(d.Σ, 1))
as(d::MvNormal{(:μ, :Λ),<:Tuple{T,C}}) where {T,C<:Cholesky} = as(Array, size(d.Λ, 1))

function as(d::MvNormal{(:σ,),Tuple{M}}) where {M<:Triangular}
    σ = d.σ
    if @inbounds all(i -> σ[i] ≠ 0, diagind(σ))
        return as(Array, size(σ, 1))
    else
        @error "Not implemented yet"
    end
end

function as(d::MvNormal{(:λ,),Tuple{M}}) where {M<:Triangular}
    λ = d.λ
    if @inbounds all(i -> λ[i] > 0, diagind(λ))
        return as(Array, size(λ, 1))
    else
        @error "Not implemented yet"
    end
end

as(::Normal) = asℝ
as(::SnedecorF) = asℝ₊
as(::StudentT) = asℝ
as(::Uniform{()}) = as𝕀
