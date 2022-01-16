
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
# Uniform over an interval
#
# We use an affine parameterization (shift, scale) versus a more traditional (lower, upper)
# because of how `asparams` works. Since `asparams` works on a per-parameter basis,
# we can enforce the `lower < upper` constraint by making `scale` be non-negative

@kwstruct Uniform(shift, scale)

function basemeasure(d::Uniform{(:shift, :scale)})
    # Lebesgue measure on (shift, scale+shift)
    lb = d.shift
    ub = d.shift + d.scale
    inbounds(x) = lb < x < ub
    constℓ = 0.0
    varℓ() = 0.0
    base = Lebesgue(ℝ)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

asparams(::Type{<:Uniform}, ::Val{:shift}) = asℝ
asparams(::Type{<:Uniform}, ::Val{:scale}) = asℝ₊

"""
The uniform probability measure on the interval (a, b).
"""
function Uniform(a, b)
    @assert(a < b)
    Uniform(shift=a, scale=b-a)
end

logdensity(d::Uniform{(:shift, :scale)}, x) = begin
    @assert (d.scale > 0.0)
    ℓ = -log(d.scale)
end

distproxy(d::Uniform{(:shift, :scale)}) = Dists.Uniform(d.shift, d.shift + d.scale)

TV.as(d::Uniform{(:shift, :scale), T}) where T = TV.ScaledShiftedLogistic(Real, d.scale, d.shift)

Base.rand(rng::AbstractRNG, T::Type, μ::Uniform{(:shift, :scale)}) = rand(rng, T)*μ.scale + μ.shift
