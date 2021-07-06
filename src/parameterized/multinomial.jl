# Multinomial distribution

using SpecialFunctions
import StatsFuns
export Multinomial


@parameterized Multinomial(n, p)

basemeasure(d::Multinomial{(:n, :p)}) = CountingMeasure(ℤ) ^ length(d.p)


@kwstruct Multinomial(n,p)

function logdensity(d::Multinomial{(:n, :p)}, x)
    p = d.p
    s = 0.0
    for j in eachindex(x)
        s += x[j] * log(p[j])
    end
    return s
end

Base.rand(rng::AbstractRNG, T::Type, μ::Multinomial) = rand(rng, Dists.Multinomial(μ.n, μ.p))

distproxy(d::Multinomial{(:p,)}) = Dists.Multinomial(d.n,d.p)


# Based on
# https://github.com/JuliaMath/Combinatorics.jl/blob/c2114a71ccfc93052efb9a9379e62b81b9388ef8/src/factorials.jl#L99
function logmultinomial(k)
    s = 0
    result = 1
    @inbounds for i in k
        s += i
        (Δresult, _) = logabsbinomial(s, i)
        result += Δresult
    end
    result
end
