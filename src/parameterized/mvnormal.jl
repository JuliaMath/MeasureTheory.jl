
# Multivariate Normal distribution

export MvNormal

function MvNormal(nt::NamedTuple{(:μ,)})
    dim = size(nt.μ)
    Affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:σ,)})
    dim = size(nt.σ, 2)
    Affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:ω,)})
    dim = size(nt.σ, 1)
    Affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:μ, :σ,)})
    dim = size(nt.σ, 2)
    Affine(nt, Normal() ^ dim)
end

function MvNormal(nt::NamedTuple{(:μ, :ω,)})
    dim = size(nt.σ, 1)
    Affine(nt, Normal() ^ dim)
end
