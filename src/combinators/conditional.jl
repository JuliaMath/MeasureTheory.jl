using NestedTuples: lazymerge

struct ConditionalMeasure{M,C} <: AbstractMeasure
    parent::M 
    constraint::C
end

Base.:|(μ::AbstractMeasure, constraint) = ConditionalMeasure(μ, constraint)

@inline function logdensity_def(cm::ConditionalMeasure, x)
    logdensity_def(cm.parent, lazymerge(cm.constraint, x))
end

@inline basemeasure(cm::ConditionalMeasure) = basemeasure(cm.parent) | cm.constraint