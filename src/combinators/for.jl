
export For
using Random
import Base

struct For{T, F, I} <: AbstractProductMeasure
    f::F
    inds::I

    @inline function For{T}(f::F, inds::I) where {T,F,I<:Tuple}
        new{T,Core.Typeof(f),I}(f, inds)
    end

    @inline For{T,F,I}(f::F, inds::I) where {T,F,I} = new{T,F,I}(f,inds)
end

@generated function For(f::F, inds::I) where {F,I<:Tuple}
    eltypes = Tuple{eltype.(I.types)...}
    quote
        $(Expr(:meta, :inline))
        T = Core.Compiler.return_type(f, $eltypes)
        For{T,F,I}(f, inds)
    end
end

TV.as(d::For) = as(Array, as(first(marginals(d))), size(first(d.inds))...)

# For(f, gen::Base.Generator) = ProductMeasure(Base.Generator(f ∘ gen.f, gen.iter))

# @inline function tailcall_iter(f, iter)
#     iter_result::Tuple = iterate(iter)
#     (val, state) = iter_result
#     tailcall_iter(f, iter, (val, state))
# end

# @inline function tailcall_iter(f, iter, ::Nothing)
#     return unit(f)
# end

# unit(+) = false

# @inline function tailcall_iter(f, iter, (val, state))
#     tailcall_iter(f, iter, (val, state), unit(f))
# end
  
# @inline function tailcall_iter(f, iter, (val, state), acc)
#     tailcall_iter(f, iter, iterate(iter, state), f(val, acc))
# end

# @inline function tailcall_iter(f, iter, ::Nothing, acc)
#     return acc
# end

# @inline function logdensity_def(d::For{T,F,I}, x::AbstractVector{X}; exec=SequentialEx(simd = true)) where {X,T,F,I<:Tuple{<:AbstractVector}}
#     f(j) = @inbounds logdensity_def(d.f(j), x[j])
#     f(j, ℓ) = ℓ + @inbounds logdensity_def(d.f(j), x[j])

#     tailcall_iter(f, eachindex(x))
# end

@inline function logdensity_def(d::For{T,F,I}, x::AbstractVector{X}; exec=SequentialEx(simd = true)) where {X,T,F,I<:Tuple{<:AbstractVector}}
    ind = only(d.inds)
    js =  eachindex(x)
    ℓ = 0
    @simd for j in js
        i = getindex(ind, j)
        Δℓ = @inbounds logdensity_def(d.f(i), x[j])
        ℓ += Δℓ
    end
    ℓ
end

function logdensity_def(d::For, x::AbstractVector; exec=SequentialEx(simd = true)) 
    @floop exec for j in eachindex(x)
        local i = (getindex(ind, j) for ind in d.inds)
        local Δℓ = @inbounds logdensity_def(d.f(i...), x[j])
        @reduce ℓ += Δℓ
    end
 
    ℓ
end
   

function logdensity_def(d::For{T,F,I}, x::AbstractArray{X}; exec=SequentialEx(simd = true)) where {T,F,I,X}
    @floop exec for j in CartesianIndices(x)
        local i = (getindex(ind, j) for ind in d.inds)
        local Δℓ = @inbounds logdensity_def(d.f(i...), x[j])
        @reduce ℓ += Δℓ
    end
 
    ℓ
end

# function logdensity_def(d::For{T,F,I}, x) where {N,T,F,I<:NTuple{N,<:Base.Generator}}
#     sum(zip(x, d.inds...)) do (xⱼ, dⱼ...)
#         logdensity_def(d.f(dⱼ...), xⱼ)
#     end
# end

# function logdensity_def(d::For{T,F,I}, x::AbstractVector) where {N,T,F,I<:NTuple{N,<:Base.Generator}}
    
#     sum(zip(x, d.inds...)) do (xⱼ, dⱼ...)
#         logdensity_def(d.f(dⱼ...), xⱼ)
#     end
# end

function marginals(d::For{T,F,Tuple{I}}) where {T,F,I}
    f = d.f
    data = first(d.inds)
    MappedArrays.ReadonlyMappedArray{T,ndims(data),typeof(data),typeof(f)}(f, data)
end

function marginals(d::For{T,F,I}) where {T,F,I}
    f = d.f
    data = d.inds
    MappedArrays.ReadonlyMultiMappedArray{T,ndims(first(data)),typeof(data),typeof(f)}(f, data)
end

function marginals(d::For{T,F,I}) where {N,T,F,I<:NTuple{N,<:Base.Generator}}
    Iterators.map(d.f, d.inds...)
end

function basemeasure(d::For{T,F,I}) where {T,F,I}
    B = typeof(basemeasure(d.f(map(first, d.inds)...)))
    sing = static(Base.issingletontype(B))
    _basemeasure(d, B, sing)
end

@inline function _basemeasure(d::For{T,F,I}, ::Type{<:WeightedMeasure{StaticFloat64{N}, B}}, ::True) where {T,F,I,N,B}
    dim = size(first(d.inds))
    weightedmeasure(static(N) * prod(dim), instance(B) ^ dim)
end


@inline function _basemeasure(d::For{T,F,I}, ::Type{<:WeightedMeasure{StaticFloat64{N}, B}}, ::False) where {T,F,I,N,B}
    dim = size(first(d.inds))
    base = For(d.inds) do j
        basemeasure(d.f(j)).base
    end
    weightedmeasure(static(N) * prod(dim), base)
end

@inline function _basemeasure(d::For{T,F,I}, ::Type{B}, ::True) where {T,F,I,B}
    instance(B) ^ size(first(d.inds))
end

@inline function _basemeasure(d::For{T,F,I}, ::Type{B}, ::False) where {T,F,I,B<:AbstractMeasure}
    f = basekleisli(d.f)
    For{B}(f, d.inds)
end

@inline function _basemeasure(d::For{T,F,I}, ::Type{B}, ::False) where {T,F,I,B}
    MeasureBase.productmeasure(basemeasure.(marginals(d)))
end

function _basemeasure(d::For{T,F,I}, ::Type{B}, ::True) where {N,T<:AbstractMeasure,F,I<:NTuple{N,<:Base.Generator},B}
    return instance(B) ^ minimum(length, d.inds)
end

function _basemeasure(d::For{T,F,I}, ::Type{B}, ::False) where {N,T<:AbstractMeasure,F,I<:NTuple{N,<:Base.Generator},B}
    f = basekleisli(d.f)
    For{B}(f, d.inds)
end

function Pretty.tile(d::For{T}) where {T}
    result = Pretty.literal("For{")
    result *= Pretty.tile(T)
    result *= Pretty.literal("}")
    result *= Pretty.list_layout(
        [
            Pretty.literal(func_string(d.f, Tuple{_eltype.(d.inds)...})),
            Pretty.tile.(d.inds)...
        ]
    )
end



"""
    For(f, base...)

`For` provides a convenient way to construct a `ProductMeasure`. There are
several options for the `base`. With Julia's `do` notation, this can look very
similar to a standard `for` loop, while maintaining semantics structure that's
easier to work with.

------------

# `For(f, base::Int...)`

When one or several `Int` values are passed for `base`, the result is treated as
depending on `CartesianIndices(base)`. 

```
julia> For(3) do λ Exponential(λ) end |> marginals
3-element mappedarray(MeasureBase.var"#17#18"{var"#15#16"}(var"#15#16"()), ::CartesianIndices{1, Tuple{Base.OneTo{Int64}}}) with eltype Exponential{(:λ,), Tuple{Int64}}:
 Exponential(λ = 1,)
 Exponential(λ = 2,)
 Exponential(λ = 3,)
```

```
julia> For(4,3) do μ,σ Normal(μ,σ) end |> marginals
4×3 mappedarray(MeasureBase.var"#17#18"{var"#11#12"}(var"#11#12"()), ::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}) with eltype Normal{(:μ, :σ), Tuple{Int64, Int64}}:
 Normal(μ = 1, σ = 1)  Normal(μ = 1, σ = 2)  Normal(μ = 1, σ = 3)
 Normal(μ = 2, σ = 1)  Normal(μ = 2, σ = 2)  Normal(μ = 2, σ = 3)
 Normal(μ = 3, σ = 1)  Normal(μ = 3, σ = 2)  Normal(μ = 3, σ = 3)
 Normal(μ = 4, σ = 1)  Normal(μ = 4, σ = 2)  Normal(μ = 4, σ = 3)
```

-------

# `For(f, base::AbstractArray...)``

In this case, `base` behaves as if the arrays are `zip`ped together before
applying the map.

```
julia> For(randn(3)) do x Exponential(x) end |> marginals
3-element mappedarray(x->Main.Exponential(x), ::Vector{Float64}) with eltype Exponential{(:λ,), Tuple{Float64}}:
 Exponential(λ = -0.268256,)
 Exponential(λ = 1.53044,)
 Exponential(λ = -1.08839,)
```

```
julia> For(1:3, 1:3) do μ,σ Normal(μ,σ) end |> marginals
3-element mappedarray((:μ, :σ)->Main.Normal(μ, σ), ::UnitRange{Int64}, ::UnitRange{Int64}) with eltype Normal{(:μ, :σ), Tuple{Int64, Int64}}:
 Normal(μ = 1, σ = 1)
 Normal(μ = 2, σ = 2)
 Normal(μ = 3, σ = 3)
```

----

# `For(f, base::Base.Generator)`

For `Generator`s, the function maps over the values of the generator:

```
julia> For(eachrow(rand(4,2))) do x Normal(x[1], x[2]) end |> marginals |> collect
4-element Vector{Normal{(:μ, :σ), Tuple{Float64, Float64}}}:
 Normal(μ = 0.255024, σ = 0.570142)
 Normal(μ = 0.970706, σ = 0.0776745)
 Normal(μ = 0.731491, σ = 0.505837)
 Normal(μ = 0.563112, σ = 0.98307)
```

"""
@inline For{T}(f, inds...) where {T} = For{T}(f, inds)
@inline For{T}(f, n::Integer) where {T} = For{T}(f, static(1):n)
@inline For{T}(f, inds::Integer...) where {T} = For{T}(i -> f(Tuple(i)...), Base.CartesianIndices(inds))

@inline For(f, inds...) = For(f, inds)
@inline For(f, n::Integer) = For(f, static(1):n)
@inline For(f, inds::Integer...) = For(i -> f(Tuple(i)...), Base.CartesianIndices(inds))
# For(f, inds::Base.Generator) = productmeasure(mymap(f, inds))

function Random.rand!(rng::AbstractRNG, d::For{T,F,I}, x) where {T,F,I} 
    mar = marginals(d)
    @inbounds for (j, dⱼ) in enumerate(mar)
        x[j] = rand(rng,dⱼ)
    end
    return x
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::For{M,F,I}) where {T,M,F,I}
    MeasureBase._rand_product(rng, T, marginals(d), M)
end