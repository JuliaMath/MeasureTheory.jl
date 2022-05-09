export ExponentialFamily

@concrete terse struct ExponentialFamily <: AbstractKleisli
    base
    mdim
    pdim
    t
    x
    a
end

function ExponentialFamily(base, mdim, pdim, t, a)
    return ExponentialFamily(base, mdim, pdim, t, I, a)
end

function MeasureBase.powermeasure(fam::ExponentialFamily, dims::NTuple{N,I}) where {N,I}
    t = Tuple((y -> f.(y) for f in fam.t))
    a(η) = BroadcastArray(fam.a, η)
    p = prod(dims)
    ExponentialFamily(fam.base^dims, fam.mdim * p, fam.pdim * p, t, fam.x, a)
end

@concrete terse struct ExpFamMeasure <: AbstractMeasure
    fam
    η # instantiated to a value
    a # instantiated to a value
end

@inline function (fam::ExponentialFamily)(β)
    η = fam.x * β
    a = fam.a(η)
    ExpFamMeasure(fam, η, a)
end

MeasureBase.basemeasure(d::ExpFamMeasure) = d.fam.base

tracedot(a::AbstractVector, b::AbstractVector) = dot(a, b)

tracedot(a::AbstractVector, x, b::AbstractVector) = dot(a, x, b)

tracedot(a, b) = sum((dot(view(a, :, j), view(b, :, j)) for j in 1:size(a, 2)))

tracedot(a, x, b) =
    sum(1:size(a, 2)) do j
        dot(view(a, :, j), x, view(b, :, j))
    end

# @inline function tracedot(a::BlockDiag, b::BlockDiag)
#     numblocks = length(a.blocks)
#     sum(tracedot(a.blocks[j], b.blocks[j]) for j in 1:length(a.blocks))
# end

# @inline function tracedot(a::BlockDiag, x::BlockDiag, b::BlockDiag)
#     numblocks = length(x.blocks)
#     sum(tracedot(a.blocks[j], x.blocks[j], b.blocks[j]) for j in 1:length(x.blocks))
# end

function logdensity_def(d::ExpFamMeasure, y)
    t = ApplyArray(vcat, (f.(y) for f in d.fam.t)...)
    η = d.η
    dot(t, η)
end

function withX(fam::ExponentialFamily, x)
    @inline t(y) = fam.t.(y)
    newx = ApplyArray(kron, x, fam.x)
    η(β) = fam.η.(β)
    a(β) = sum(fam.a, β)
    ExponentialFamily(fam.base^size(x, 1), t, x, η, a)
end

@concrete terse struct ExpFamLikelihood <: AbstractLikelihood
    fam
    y
    tᵀx
    c
end

export likelihood

function regression(fam, uᵀ, vᵀ)
end

function likelihood(fam::ExponentialFamily, y)
    c = logdensityof(fam.base, y)
    t = ApplyArray(vcat, (f.(y) for f in fam.t)...)
    tᵀx = t' * fam.x
    ExpFamLikelihood(fam, y, tᵀx, c)
end

@inline function logdensity_def(ℓ::ExpFamLikelihood, β)
    xβ = ApplyArray(*, ℓ.fam.x, β)
    a = sum(ℓ.fam.a(xβ))
    # a = sum(ℓ.fam.a, ApplyArray(*, ℓ.fam.uᵀ', ℓ.fam.vᵀ, β))
    ℓ.c + dot(ℓ.tᵀx, β) - a
end

basemeasure(fam::ExponentialFamily) = fam.base

# function stack_functions(funs, inds)
#     function(x::AbstractArray{T,N}) where {T,N}
#         ApplyArray(cat, )
# end
