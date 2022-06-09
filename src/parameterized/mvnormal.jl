
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

@useproxy MvNormal

proxy(d::MvNormal) = affine(params(d), Normal()^supportdim(d))

rand(rng::AbstractRNG, ::Type{T}, d::MvNormal) where {T} = rand(rng, T, proxy(d))

insupport(d::MvNormal, x) = insupport(proxy(d), x)

proxy(d::MvNormal{(:Σ,),Tuple{C}}) where {C<:Cholesky}    = MvNormal(σ = getL(d.Σ))
proxy(d::MvNormal{(:Λ,),Tuple{C}}) where {C<:Cholesky}    = MvNormal(λ = getL(d.λ))
proxy(d::MvNormal{(:μ, :Σ),Tuple{C}}) where {C<:Cholesky} = MvNormal(μ = d.μ, σ = getL(d.Σ))
proxy(d::MvNormal{(:μ, :Λ),Tuple{C}}) where {C<:Cholesky} = MvNormal(μ = d.μ, λ = getL(d.Λ))

# function MvNormal(nt::NamedTuple{(:μ,)})
#     dim = size(nt.μ)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:σ,)})
#     dim = colsize(nt.σ)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:λ,)})
#     dim = rowsize(nt.λ)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:μ, :σ,)})
#     dim = colsize(nt.σ)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:μ, :λ,)})
#     dim = rowsize(nt.λ)
#     affine(nt, Normal() ^ dim)
# end
