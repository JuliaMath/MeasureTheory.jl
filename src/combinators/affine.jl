asparams(::Affine, ::Val{:μ}) = asℝ
asparams(::Affine, ::Val{:σ}) = asℝ₊
asparams(::Type{A}, ::Val{:μ}) where {A<:Affine} = asℝ
asparams(::Type{A}, ::Val{:σ}) where {A<:Affine} = asℝ₊

asparams(::Affine{N,M,T}, ::Val{:μ}) where {N,M,T<:AbstractArray} =
    as(Array, asℝ, size(d.μ))
asparams(::Affine{N,M,T}, ::Val{:σ}) where {N,M,T<:AbstractArray} =
    as(Array, asℝ, size(d.σ))
