
# Multivariate Normal distribution

import StatsFuns
using LinearAlgebra
export MvNormal
using Random
import Base
using Tullio
using LoopVectorization

@parameterized MvNormal(μ,Σ)

@kwstruct MvNormal(Σ)

@kwstruct MvNormal(μ,Σ)

@kwstruct MvNormal(L)

@kwstruct MvNormal(μ,L)

@kwstruct MvNormal(μ,Ω)

@kwstruct MvNormal(Ω)


function logdensity(d::MvNormal{(:L2,)}, y::AbstractVector{T}) where {T}
    L = d.L2
    k = first(size(L))
    @assert length(y) == k

    # if z = Lᵗy, the logdensity depends on `det(L)` and `z ⋅ z`. So we find `z`
    z_dot_z = zero(T)
    @simd for j ∈ 1:k
        zj = zero(T)
        for i ∈ j:k
            @inbounds zj += L[i,j] * y[i]
        end
        z_dot_z += zj ^ 2
    end

    -k/2 * log2π + logdet_pos(L) - z_dot_z/2
end

@kwstruct MvNormal(L2)

using StrideArrays
using StaticArrays


@inline function logdet_pos(A::Union{UpperTriangular{T},LowerTriangular{T}}) where T
  abs_det = zero(real(T))
  @turbo for i in 1:size(A,1)
    diag_i = A.data[i,i]
    abs_det += log(diag_i)
  end
  return abs_det
end

function logdensity(d::MvNormal{(:μ,:L)}, y::AbstractArray{T}) where {T}
    x = StrideArray{T}(undef, size(y))
    @inbounds for j in eachindex(y)
        x[j] = y[j] - d.μ[j]
    end
    GC.@preserve x logdensity(MvNormal(L=d.L), x)
end

@inline function logdensity(d::MvNormal{(:L,)}, y::AbstractVector{T}) where {T}
    L = d.L
    k = length(y)
    z = StrideArray{T}(undef, (k,))

    # Solve `y = Lz` for `z`. We need this only as a way to calculate `z ⋅ z`
    z_dot_z = zero(T)
    @inbounds for i ∈ 1:k
        tmp = zero(T)
        for j in 1:(i-1)
            tmp += L[i,j] * z[j]
        end
        zi = (y[i] - tmp) / L[i,i]
        z_dot_z += zi ^ 2
        z[i] = zi
    end

    -k/2 * log2π - logdet_pos(L) - z_dot_z/2
end

@generated function logdensity(d::MvNormal{(:L,)}, y::SizedVector{K,T}) where {K,T}
  log2π = log(big(2)*π)
  if K > 16
    quote
      L = parent(d.L)
      k = StaticInt($K)

      z = StrideArray{$T}(undef, (k,))

      # Solve `y = Lz` for `z`. We need this only as a way to calculate `z ⋅ z`
      z_dot_z = zero($T)
      @inbounds @fastmath for i ∈ 1:k
        tmp = zero($T)
        for j in 1:(i-1)
          tmp += L[i,j] * z[j]
        end
        zi = (y[i] - tmp) / L[i,i]
        z_dot_z += zi ^ 2
        z[i] = zi
      end

      $(T(-K/2 * log2π)) - logdet_pos(L) - z_dot_z/2
    end
  else # replace `for i in 1:K` with `Base.Cartesian.@nexprs $K i -> begin`
    quote
      L = d.L
      # Solve `y = Lz` for `z`. We need this only as a way to calculate `z ⋅ z`
      z_dot_z = zero($T)

      @inbounds begin # `@fastmath` needs to be moved inside the `@nexprs`
        Base.Cartesian.@nexprs $K i -> begin
          tmp_i = zero($T)
          Base.Cartesian.@nexprs i-1 j -> begin
            @fastmath tmp_i += L[i,j] * z_j
          end
          @fastmath z_i = (y[i] - tmp_i) / L[i,i]
          @fastmath z_dot_z += z_i ^ 2
        end
      end

      $(T(-K/2 * log2π)) - logdet_pos(L) - z_dot_z/2
    end
  end
end

# S = @MMatrix(randn(10,15)) |> x -> Symmetric(x * x',:L);
# d = (L = cholesky(S).L,)
# y = SizedVector{10}(randn(10));

# z = StrideArray{Float64}(undef, (StaticInt(3),))
# b = z.data

# function MeasureTheory.basemeasure(μ::MvNormal{N, T, I,I}) where {N, T, I}
#     return (1 / sqrt2π)^I * Lebesgue(ℝ)^I
# end

# sampletype(d::MvNormal{N, T, I, J}) where {N,T,I,J} = Vector{Float64}



# function Random.rand!(rng::AbstractRNG, d::MvNormal{(:A, :b),T,I,J}, x::AbstractArray) where {T,I,J}
#     z = getcache(d, J)
#     rand!(rng, Normal()^J, z)
#     copyto!(x, d.b)
#     mul!(x, d.A, z, 1.0, 1.0)
#     return x
# end

# function logdensity(d::MvNormal{(:A,:b)}, x)
#     cache = getcache(d, size(d))
#     z = d.A \ (x - d.b)
#     return logdensity(MvNormal(), z) - logabsdet(d.A)
# end

# ≪(::MvNormal, ::Lebesgue{ℝ}) = true

# function logdensity(d::MvNormal{(:Σ⁻¹,)}, x)
#     @tullio ℓ = -0.5 * x[i] * d.Σ⁻¹[i,j] * x[j]
#     return ℓ
# end

# mvnormaldims(nt::NamedTuple{(:Σ⁻¹,)}) = size(nt.Σ⁻¹)
