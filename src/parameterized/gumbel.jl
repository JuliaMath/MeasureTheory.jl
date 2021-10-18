# Gumbel distribution

export Gumbel

@parameterized Gumbel() 

basemeasure(::Gumbel{()}) = Lebesgue(ℝ)

@kwstruct Gumbel()

@kwstruct Gumbel(β)

@kwstruct Gumbel(μ,β)

# map affine names to those more common for Gumbel
for N in [(:μ,), (:σ,), (:μ,:σ)]
    G = replace(N, :σ => :β)
    @eval begin
        proxy(d::Gumbel{$G}) = affine(NamedTuple{$N}(values(params(d))), Gumbel())
        logdensity(d::Gumbel{$G}, x) = logdensity(proxy(d), x)
        basemeasure(d::Gumbel{$G}) = basemeasure(proxy(d))
    end
end

function logdensity(d::Gumbel{()} , x)
    return -exp(-x) - x
end

import Base

function Base.rand(rng::AbstractRNG, d::Gumbel{()})
    u = rand(rng)
    -log(-log(u))
end

TV.as(::Gumbel) = asℝ

≪(::Gumbel, ::Lebesgue{X}) where X <: Real = true

distproxy(::Gumbel{()}) = Dists.Gumbel()
