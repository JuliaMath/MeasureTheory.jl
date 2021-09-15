asparams(::Affine, ::Val{:μ}) = asℝ
asparams(::Affine, ::Val{:σ}) = asℝ₊
asparams(::Type{A}, ::Val{:μ}) where {A<:Affine} = asℝ
asparams(::Type{A}, ::Val{:σ}) where {A<:Affine} = asℝ₊

