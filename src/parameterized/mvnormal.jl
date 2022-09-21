
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



for N in setdiff(AFFINEPARS, [(:μ,)])
    @eval begin
        function as(d::MvNormal{$N})
            p = proxy(d)
            if rank(getfield(p, :f)) == only(supportdim(d))
                return as(Array, supportdim(d))
            else
                @error "Not yet implemented"
            end
        end
    end
end

supportdim(d::MvNormal) = supportdim(params(d))

supportdim(nt::NamedTuple{(:Σ,)}) = size(nt.Σ, 1)
supportdim(nt::NamedTuple{(:μ, :Σ)}) = size(nt.Σ, 1)
supportdim(nt::NamedTuple{(:Λ,)}) = size(nt.Λ, 1)
supportdim(nt::NamedTuple{(:μ, :Λ)}) = size(nt.Λ, 1)

@useproxy MvNormal

for N in [(:Σ,), (:μ, :Σ), (:Λ,), (:μ, :Λ)]
    @eval basemeasure_depth(d::MvNormal{$N}) = static(2)
end

proxy(d::MvNormal) = affine(params(d), Normal()^supportdim(d))

rand(rng::AbstractRNG, ::Type{T}, d::MvNormal) where {T} = rand(rng, T, proxy(d))

insupport(d::MvNormal, x) = insupport(proxy(d), x)

# Note: (C::Cholesky).L may or may not make a copy, depending on C.uplo, which is not included in the type
@inline function proxy(d::MvNormal{(:Σ,),Tuple{C}}) where {C<:Cholesky}
    affine((σ = d.Σ.L,), Normal()^supportdim(d))
end
@inline function proxy(d::MvNormal{(:Λ,),Tuple{C}}) where {C<:Cholesky}
    affine((λ = d.Λ.L,), Normal()^supportdim(d))
end
@inline function proxy(d::MvNormal{(:μ, :Σ),Tuple{T,C}}) where {T,C<:Cholesky}
    affine((μ = d.μ, σ = d.Σ.L), Normal()^supportdim(d))
end
@inline function proxy(d::MvNormal{(:μ, :Λ),Tuple{T,C}}) where {T,C<:Cholesky}
    affine((μ = d.μ, λ = d.Λ.L), Normal()^supportdim(d))
end
