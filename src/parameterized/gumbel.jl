# Gumbel distribution

export Gumbel

@parameterized Gumbel()

basemeasure(::Gumbel{()}) = LebesgueBase()

@kwstruct Gumbel()

@kwstruct Gumbel(β)

@kwstruct Gumbel(μ, β)

# map affine names to those more common for Gumbel
for N in [(:μ,), (:σ,), (:μ, :σ)]
    G = tuple(replace(collect(N), :σ => :β)...)
    @eval begin
        proxy(d::Gumbel{$G}) = affine(NamedTuple{$N}(values(params(d))), Gumbel())
        logdensity_def(d::Gumbel{$G}, x) = logdensity_def(proxy(d), x)
        basemeasure(d::Gumbel{$G}) = basemeasure(proxy(d))
    end
end

@inline function logdensity_def(d::Gumbel{()}, x)
    return -exp(-x) - x
end

import Base

function Base.rand(rng::AbstractRNG, ::Type{T}, d::Gumbel{()}) where {T}
    u = rand(rng, T)
    -log(-log(u))
end

≪(::Gumbel, ::Lebesgue{X}) where {X<:Real} = true

insupport(::Gumbel, x) = true

proxy(::Gumbel{()}) = Dists.Gumbel()
