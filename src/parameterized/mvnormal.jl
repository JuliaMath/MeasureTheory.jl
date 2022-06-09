
# Multivariate Normal distribution

export MvNormal

@parameterized MvNormal(μ, σ)

# MvNormal(;kwargs...) = MvNormal(kwargs)

@kwstruct MvNormal(μ)

@kwstruct MvNormal(σ)
@kwstruct MvNormal(λ)

@kwstruct MvNormal(μ, σ)
@kwstruct MvNormal(μ, λ)

@kwstruct MvNormal(Σ)
@kwstruct MvNormal(Λ)
@kwstruct MvNormal(μ, Σ)
@kwstruct MvNormal(μ, Λ)

as(d::MvNormal{(:μ,)}) = as(Array, length(d.μ))

as(d::MvNormal{(:Σ,),Tuple{C}}) where {C<:Cholesky}    = as(Array, size(d.Σ, 1))
as(d::MvNormal{(:Λ,),Tuple{C}}) where {C<:Cholesky}    = as(Array, size(d.Λ, 1))
as(d::MvNormal{(:μ, :Σ),Tuple{C}}) where {C<:Cholesky} = as(Array, size(d.Σ, 1))
as(d::MvNormal{(:μ, :Λ),Tuple{C}}) where {C<:Cholesky} = as(Array, size(d.Λ, 1))

function as(d::MvNormal{(:σ,),Tuple{M}}) where {M<:Triangular}
    σ = d.σ
    if @inbounds all(i -> σ[i] > 0, diagind(σ))
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

for N in setdiff(AFFINEPARS, [(:μ,)])
    @eval begin
        function as(d::MvNormal{$N})
            p = proxy(d)
            if rank(getfield(p,:f)) == only(supportdim(d))
                return as(Array, supportdim(d))
            else
                @error "Not yet implemented"
            end
        end
    end
end

supportdim(d::MvNormal) = supportdim(params(d))

supportdim(nt::NamedTuple{(:Σ,)}) = size(nt.Σ, 1)
supportdim(nt::NamedTuple{(:μ,:Σ)}) = size(nt.Σ, 1)
supportdim(nt::NamedTuple{(:Λ,)}) = size(nt.Λ, 1)
supportdim(nt::NamedTuple{(:μ,:Λ)}) = size(nt.Λ, 1)

@useproxy MvNormal

for N in [(:Σ,), (:μ,:Σ), (:Λ,), (:μ,:Λ)]
    @eval basemeasure_depth(d::MvNormal{$N}) = static(2)
end

proxy(d::MvNormal) = affine(params(d), Normal()^supportdim(d))

rand(rng::AbstractRNG, ::Type{T}, d::MvNormal) where {T} = rand(rng, T, proxy(d))

insupport(d::MvNormal, x) = insupport(proxy(d), x)

# Note: (C::Cholesky).L may or may not make a copy, depending on C.uplo, which is not included in the type
@inline proxy(d::MvNormal{(:Σ,),Tuple{C}}) where {C<:Cholesky}    = affine((σ = d.Σ.L,), Normal()^supportdim(d))
@inline proxy(d::MvNormal{(:Λ,),Tuple{C}}) where {C<:Cholesky}    = affine((λ = d.Λ.L,), Normal()^supportdim(d))
@inline proxy(d::MvNormal{(:μ, :Σ),Tuple{C}}) where {C<:Cholesky} = affine((μ = d.μ, σ = d.Σ.L), Normal()^supportdim(d))
@inline proxy(d::MvNormal{(:μ, :Λ),Tuple{C}}) where {C<:Cholesky} = affine((μ = d.μ, λ = d.Λ.L), Normal()^supportdim(d))
