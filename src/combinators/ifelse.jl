struct IfElseMeasure{P,T,F}
    p::P
    t::T
    f::F
end

insupport(d::IfElseMeasure, x) = insupport(d.t, x) || insupport(d.f, x)

function logdensity_def(d::IfElseMeasure, x)
    logdensity_def(d.p * d.t + (1 - d.p) * d.f, x)
end

basemeasure(d::IfElseMeasure) = d.t + d.f

IfElse.ifelse(p::Bernoulli, t, f) = IfElseMeasure(mean(p), t, f)