
# Multivariate Normal distribution

export MvNormal

function MvNormal(nt::NamedTuple{(:μ,)})
    dim = size(nt.μ)
    affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:σ,)})
    dim = colsize(nt.σ)
    affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:ω,)})
    dim = rowsize(nt.ω)
    affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:μ, :σ,)})
    dim = colsize(nt.σ)
    affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:μ, :ω,)})
    dim = rowsize(nt.ω)
    affine(nt, Normal() ^ dim)
end
