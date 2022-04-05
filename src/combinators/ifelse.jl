struct IfElseMeasure{B,T,F} <: AbstractMeasure
    b::B
    t::T
    f::F
end

insupport(d::IfElseMeasure, x) = insupport(d.t, x) || insupport(d.f, x)

function logdensity_def(d::IfElseMeasure, x)
    p = mean(d.b)
    logdensity_def(p * d.t + (1 - p) * d.f, x)
end

basemeasure(d::IfElseMeasure) = d.t + d.f

@inline function Base.rand(rng::AbstractRNG, ::Type{T}, m::IfElseMeasure) where {T}
    c = rand(rng, T, m.b)
    if c
        return rand(rng, T, m.t)
    else
        return rand(rng, T, m.f)
    end
end


IfElse.ifelse(b::Bernoulli, t, f) = IfElseMeasure(b, t, f)