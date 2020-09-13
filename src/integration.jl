import Distributions

const Dists=Distributions

struct Integral{F,M}
    f::F
    μ::M
end

∫(μ) = ∫(identity,μ)
∫(f, μ) = Integral(f, μ)

∫(::typeof(identity), ::Dists.Distribution) = 1.0

∫(Dists.Normal()) do x x^2 end