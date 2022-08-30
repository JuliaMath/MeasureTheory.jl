function convolve(μ::Normal, ν::Normal)
    Normal(mean(μ) + mean(ν), hypot(std(μ), std(ν)))
end

function bind(μ::Normal,
    ::ParameterizedTransitionKernel{Type{Normal{(:μ,)}}, typeof(identity), (:μ,), Tuple{typeof(identity)}})
   convolve(μ, Normal())
end

function bind(μ::Normal,
    k::ParameterizedTransitionKernel{Type{Normal{(:μ, :σ)}}, typeof(identity), (:μ, :σ), Tuple{typeof(identity), T}} where T<:Number)
   convolve(μ, Normal(σ=k.param_maps.σ))
end
