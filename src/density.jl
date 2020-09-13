struct Density{M,B}
    μ::M
    base::B
end

struct LogDensity{M,B}
    μ::M
    base::B
end


density(μ, base) = Density(μ, base)
logdensity(μ, base) = LogDensity(μ, base)

density(μ::Distribution{Univariate,Continuous}, x::Real) = pdf(μ,x)
logdensity(μ::Distribution{Univariate,Continuous}, x::Real) = logpdf(μ,x)
