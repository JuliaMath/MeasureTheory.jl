
# Uniform distribution

export Uniform

@parameterized Uniform()
@kwstruct Uniform()


###############################################################################
# Standard Uniform

function basemeasure(::Uniform{()})
    inbounds(x) = 0 < x < 1
    constℓ = 0.0
    varℓ() = 0.0
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

distproxy(::Uniform{()}) = Dists.Uniform()

logdensity(d::Uniform{()}, x) = 0.0

TV.as(::Uniform{()}) = as𝕀

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{()}) = rand(rng, T)

###############################################################################
# Uniform
#
# The parameterization was chosen to enforce the constraint that the parameters
# define an interval. Since `asparams` operates on parameters independently,
# we parameterize with the left point of the interval (lower bound) and the width
# of the interval, which is parameterized as positive real number (`asℝ₊`).

@kwstruct Uniform(left, width)

function basemeasure(d::Uniform{(:left, :width)})
    # Lebesgue measure on (left, width+left)
    lb = d.left
    ub = d.left + d.width
    inbounds(x) = lb < x < ub
    constℓ = 0.0
    varℓ() = 0.0
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

asparams(::Type{<:Uniform{(:left, :width)}}, ::Val{:left}) = asℝ
asparams(::Type{<:Uniform{(:left, :width)}}, ::Val{:width}) = asℝ₊

"""
The uniform probability measure on the interval (a, b).
"""
function Uniform(a, b)
    @assert(a < b)
    Uniform(left=a, width=b-a)
end

logdensity(d::Uniform{(:left, :width)}, x) = -log(d.width)

distproxy(d::Uniform{(:left, :width)}) = Dists.Uniform(d.left, d.left + d.width)

TV.as(d::Uniform{(:left, :width), T}) where T = TV.as(Real, d.left, d.left + d.width)

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{(:left, :width)}) = rand(rng, T)*μ.width + μ.left
