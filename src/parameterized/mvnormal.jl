
# Multivariate Normal distribution

import StatsFuns
using LinearAlgebra
export MvNormal
using Random
import Base
using Tullio
using LoopVectorization

using StrideArrays
using StaticArrays

@parameterized MvNormal(μ, Σ)

@kwstruct MvNormal(Σ)

@kwstruct MvNormal(μ, Σ)

@kwstruct MvNormal(σ)

@kwstruct MvNormal(μ, σ)

@kwstruct MvNormal(μ, Ω)

@kwstruct MvNormal(Ω)

@kwstruct MvNormal(ω)

const LowerCholesky{T} = Cholesky{T, <:LowerTriangular} 
const UpperCholesky{T} = Cholesky{T, <:UpperTriangular} 

@inline function logdet_pos(C::Cholesky)
    logdet_pos(getfield(C, :factors))
end


@inline function logdet_pos(A::Union{UpperTriangular{T},LowerTriangular{T}}) where {T}
    # ∑ᵢ log(aᵢ * 2^bᵢ) = log(∏ᵢ aᵢ) + log2 * ∑ᵢ bᵢ

    prod_ai = one(real(T))
    sum_bi = zero(real(T))

    @turbo for i = 1:ArrayInterface.size(A,1)
        diag_i = A.data[i, i]
        bi = exponent(diag_i)
        sum_bi += bi
        ai = significand(diag_i)
        prod_ai *= ai
    end
    return log(prod_ai) + logtwo * sum_bi      
end


###############################################################################
# MvNormal(σ)

logdensity(d::MvNormal{(:σ,), <:Tuple{<:Cholesky}}, y) = logdensity_mvnormal_σ(getfield(d.σ, :factors), y)

@generated function logdensity(
    d::MvNormal{(:σ,),Tuple{LowerTriangular{T2,S}}},
    y::AbstractArray{T},
) where {k,T,T2,S<:SizedMatrix{k,k}}
    log2π = log(big(2) * π)
    k = StaticInt(k)
    if k > 16
        quote
            σ = d.σ
            z = StrideArray{$T}(undef, ($k,))

            # Solve `y = σz` for `z`. We need this only as a way to calculate `z ⋅ z`
            z_dot_z = zero($T)
            @inbounds @fastmath for i ∈ 1:$k
                tmp = zero($T)
                for j = 1:(i-1)
                    tmp += σ[i, j] * z[j]
                end
                zi = (y[i] - tmp) / σ[i, i]
                z_dot_z += zi^2
                z[i] = zi
            end

            $(T(-k / 2 * log2π)) - logdet_pos(σ) - z_dot_z / 2
        end
    else # replace `for i in 1:K` with `Base.Cartesian.@nexprs $K i -> begin`
        quote
            σ = d.σ
            # Solve `y = σz` for `z`. We need this only as a way to calculate `z ⋅ z`
            z_dot_z = zero($T)

            @inbounds begin # `@fastmath` needs to be moved inside the `@nexprs`
                Base.Cartesian.@nexprs $k i -> begin
                    tmp_i = zero($T)
                    Base.Cartesian.@nexprs i - 1 j -> begin
                        @fastmath tmp_i += σ[i, j] * z_j
                    end
                    @fastmath z_i = (y[i] - tmp_i) / σ[i, i]
                    @fastmath z_dot_z += z_i^2
                end
            end

            $(T(-k / 2 * log2π)) - logdet_pos(σ) - z_dot_z / 2
        end
    end
end

@generated function logdensity2(
    d::MvNormal{(:σ,),Tuple{LowerTriangular{T2,S}}},
    y::AbstractArray{T},
) where {k,T,T2,S<:StaticMatrix{k,k}}
    log2π = log(big(2) * π)
    if k > 16
        quote
            # k = StaticInt($K)
            σ = d.σ
            P = parent(σ)
            tmp = StrideArray{$T}(undef, (StaticInt{$k}(),))
            @turbo for i ∈ eachindex(tmp)
                tmp[i] = zero($T)
            end
            # Solve `y = σz` for `z`. We need this only as a way to calculate `z ⋅ z`
            z_dot_z = zero($T)
            @inbounds for i ∈ 1:$k
                @fastmath z_i = (y[i] - tmp[i]) / σ[i, i]
                @fastmath z_dot_z += z_i * z_i
                @turbo check_empty = true for j ∈ i+1:$k
                    tmp[j] += P[j, i] * z_i
                end
            end
            $(T(-k / 2 * log2π)) - logdet_pos(σ) - z_dot_z / 2
        end
    else # replace `for i in 1:K` with `Base.Cartesian.@nexprs $K i -> begin`
        quote
            σ = d.σ
            # Solve `y = σz` for `z`. We need this only as a way to calculate `z ⋅ z`
            z_dot_z = zero($T)

            @inbounds begin # `@fastmath` needs to be moved inside the `@nexprs`
                Base.Cartesian.@nexprs $k i -> tmp_i = zero($T)
                Base.Cartesian.@nexprs $k i -> begin
                    # tmp_i = zero($T)
                    # Base.Cartesian.@nexprs i - 1 j -> begin
                    # @fastmath tmp_i += σ[i, j] * z_j
                    # end
                    @fastmath z_i = (y[i] - tmp_i) / σ[i, i]
                    Base.Cartesian.@nexprs $k - i j -> begin
                        @fastmath tmp_{j + i} += σ[j+i, i] * z_i
                    end
                    @fastmath z_dot_z += z_i^2
                end
            end

            $(T(-k / 2 * log2π)) - logdet_pos(σ) - z_dot_z / 2
        end
    end
end


###############################################################################
# MvNormal(ω)



logdensity(d::MvNormal{(:ω,), <:Tuple{<:Cholesky}}, y) = logdensity_mvnormal_ω(getfield(d.ω, :factors), y)

logdensity_mvnormal_ω(UL::Union{UpperTriangular, LowerTriangular}, y::AbstractVector) = logdensity_mvnormal_ω(KnownSize(UL), y)

@inline @generated function logdensity_mvnormal_ω(U::KnownSize{Tuple{k,k}, <:UpperTriangular}, y::AbstractVector{T}) where {k,T}
    log2π = log(big(2) * π)
    quote
        U = U.value
        Udata = U.data
        # if z = Lᵗy, the logdensity depends on `det(U)` and `z ⋅ z`. So we find `z`
        z_dot_z = zero(T)
        @fastmath for j ∈ 1:$k
            zj = zero(T)
            for i ∈ 1:j
                @inbounds zj += Udata[i, j] * y[i]
            end
            z_dot_z += zj^2
        end

        $(T(-k / 2 * log2π)) + logdet_pos(U) - z_dot_z / 2
    end
end

@inline function logdensity_mvnormal_ω(U::KnownSize{Tuple{nothing, nothing}, <:UpperTriangular}, y::AbstractVector{T}) where {T}
    k = first(size(U))
    @assert length(y) == k

    # if z = Lᵗy, the logdensity depends on `det(U)` and `z ⋅ z`. So we find `z`
    z_dot_z = zero(T)
    for j ∈ 1:k
        zj = zero(T)
        for i ∈ 1:j
            @inbounds zj += U[i, j] * y[i]
        end
        z_dot_z += zj^2
    end

    -k / 2 * log2π + logdet_pos(U) - z_dot_z / 2
end






function logdensity(d::MvNormal{(:μ, :σ)}, y::AbstractArray{T}) where {T}
    x = StrideArray{T}(undef, ArrayInterface.size(y))
    @inbounds for j in eachindex(y)
        x[j] = y[j] - d.μ[j]
    end
    GC.@preserve x logdensity(MvNormal(σ = d.σ), x)
end
using StrideArrays, StaticArrays, LoopVectorization, LinearAlgebra





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

@testset "MvNormal" begin
    @testset "MvNormal(L)" begin
        n = 5
        t = CorrCholesky(5)
        L = transform(t, randn(dimension(t)))

    end
end
