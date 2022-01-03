
# Multivariate Normal distribution

export MvNormal

@parameterized MvNormal(μ, σ)

# MvNormal(;kwargs...) = MvNormal(kwargs)

@kwstruct MvNormal(μ)
@kwstruct MvNormal(σ)
@kwstruct MvNormal(ω)
@kwstruct MvNormal(μ, σ)
@kwstruct MvNormal(μ, ω)

supportdim(d::MvNormal) = supportdim(params(d))

@useproxy MvNormal

proxy(d::MvNormal) = affine(params(d), Normal()^supportdim(d))

rand(rng::AbstractRNG, ::Type{T}, d::MvNormal) where {T} = rand(rng, T, proxy(d))

# function MvNormal(nt::NamedTuple{(:μ,)})
#     dim = size(nt.μ)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:σ,)})
#     dim = colsize(nt.σ)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:ω,)})
#     dim = rowsize(nt.ω)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:μ, :σ,)})
#     dim = colsize(nt.σ)
#     affine(nt, Normal() ^ dim)
# end

# function MvNormal(nt::NamedTuple{(:μ, :ω,)})
#     dim = rowsize(nt.ω)
#     affine(nt, Normal() ^ dim)
# end
