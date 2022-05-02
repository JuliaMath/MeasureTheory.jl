# InverseGaussian distribution

export InverseGaussian

@parameterized InverseGaussian()

@kwstruct InverseGaussian()
@kwstruct InverseGaussian(μ)
@kwstruct InverseGaussian(μ,λ)

InverseGaussian(μ) = InverseGaussian((μ=μ,))

InverseGaussian(μ, λ) = InverseGaussian((μ=μ,λ=λ))

@useproxy InverseGaussian{(:k, :θ)}

function logdensity_def(d::InverseGaussian{()}, x) 
    return (-3log(x) - (x - 1)^2 / x)/2
end


function logdensity_def(d::InverseGaussian{(:μ,)}, x) 
    μ= d.μ
    return (-3log(x) - (x - μ)^2 / (μ^2 * x))/2
end


function logdensity_def(d::InverseGaussian{(:μ,:λ)}, x) 
    μ, λ = d.μ, d.λ
    return (log(λ) - 3log(x) - λ * (x - μ)^2 / (μ^2 * x))/2
end

function basemeasure(::InverseGaussian)
    ℓ = static(- log2π / 2)
    weightedmeasure(ℓ, Lebesgue())
end

Base.rand(rng::AbstractRNG, T::Type, d::InverseGaussian) =
    rand(rng, proxy(d))


TV.as(::InverseGaussian) = asℝ₊

insupport(::InverseGaussian, x) = x > 0

@kwstruct InverseGaussian(μ, ϕ)

function logdensity_def(d::InverseGaussian{(:μ,:ϕ)}, x) 
    μ, ϕ = d.μ, d.ϕ
    return ( - 3log(x) - (x - μ)^2 / (ϕ * μ^2 * x))/2
end

function basemeasure(d::InverseGaussian{(:μ,:ϕ)})
    ℓ = static(-0.5) * (static(log2π) + log(d.ϕ))
    weightedmeasure(ℓ, Lebesgue())
end