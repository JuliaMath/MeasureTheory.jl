# export Mixture

using LinearAlgebra
@parameterized Mixture(components, logweights)


function logdens2(d::Mixture{(:components, :logweights)}, x::T) where {T}

    components = d.components
    logweights = d.logweights

    # Ref: [https://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html](https://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html)

    r, α = zero(T), T(-Inf)
    s, β = zero(T), T(-Inf)

    @inbounds @fastmath for i in eachindex(components)
        dᵢ = components[i]
        logwᵢ = logweights[i]
        ℓᵢ = logdensity(dᵢ, x) + logwᵢ
        
        if ℓᵢ ≤ α
            r += exp(ℓᵢ - α)
        else
            r *= exp(α - ℓᵢ)
            r += one(T)
            α = ℓᵢ
        end

        if logwᵢ ≤ β
            s += exp(logwᵢ - β)
        else
            s *= exp(β - logwᵢ)
            s += one(T)
            β = logwᵢ
        end
    end

    # ℓ = log(r) + α
    # logtotal = log(s) + β

    return log(r/s) + α - β
end



# function logdens1(d::Mixture{(:components, :logweights)}, x::T) where {T}

#     components = d.components
#     logweights = d.logweights

#     # Ref: [https://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html](https://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html)

#     r, α = zero(T), T(-Inf)
#     s, β = zero(T), T(-Inf)

#     @inbounds for i in eachindex(components, logweights)
#         dᵢ = components[i]
#         logwᵢ = logweights[i]
#         ℓᵢ = logdensity(dᵢ, x) + logwᵢ
        
#         if ℓᵢ ≤ α
#             r += exp(ℓᵢ - α)
#         else
#             r *= exp(α - ℓᵢ)
#             r += one(T)
#             α = ℓᵢ
#         end

#         if logwᵢ ≤ β
#             s += exp(logwᵢ - β)
#         else
#             s *= exp(β - logwᵢ)
#             s += one(T)
#             β = logwᵢ
#         end
#     end

#     # ℓ = log(r) + α
#     # logtotal = log(s) + β

#     return log(r/s) + α - β
# end

using StatsFuns: logsumexp

function logdensity(d::Mixture{(:components, :logweights)}, x::T) where {T}
    lp = @inbounds logsumexp((lw + logdensity(d.components[i], x) for (i,lw) in enumerate(d.logweights)))
    return lp - logsumexp(d.logweights)
end

logp = [1.0, -1.0]
p = exp.(logp)

m = Mixture((Normal(), Cauchy()), tuple(logp...));
d = Dists.MixtureModel([Dists.Normal(), Dists.Cauchy()], normalize(p,1))


logdens2(m, 0.2)
logdens3(m, 0.2)

using BenchmarkTools

@btime logdens2($m, 0.2)
@btime logdens3($m, 0.2)
