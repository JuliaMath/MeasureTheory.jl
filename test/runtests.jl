using MeasureTheory
using Test
using StatsFuns
using TransformVariables: transform
using Base.Iterators: take
using Random
using LinearAlgebra
using DynamicIterators: trace, TimeLift

function draw2(μ)
    x = rand(μ)
    y = rand(μ)
    while x == y
        y = rand(μ)
    end
    return (x,y)
end

@testset "Parameterized Measures" begin
    @testset "Binomial" begin
        D = Binomial{(:n, :p)}
        par = merge((n=20,),transform(asparams(D, (n=20,)), randn(1)))
        d = D(par)
        (n,p) = (par.n, par.p)
        logitp = logit(p)
        probitp = norminvcdf(p)
        y = rand(d)

        ℓ = logdensity(Binomial(;n, p), y)
        @test ℓ ≈ logdensity(Binomial(;n, logitp), y)
        @test ℓ ≈ logdensity(Binomial(;n, probitp), y)

        @test_broken logdensity(Binomial(n,p), CountingMeasure(ℤ[0:n]), x) ≈ binomlogpdf(n,p,x)
    end

    @testset "Normal" begin
        D = Normal{(:μ,:σ)}
        par = transform(asparams(D), randn(2))
        d = D(par)
        @test params(d) == par

        μ = par.μ
        σ = par.σ
        σ² = σ^2
        τ = 1/σ²
        logσ = log(σ)
        y = rand(d)

        ℓ = logdensity(Normal(;μ,σ), y)
        @test ℓ ≈ logdensity(Normal(;μ,σ²), y)
        @test ℓ ≈ logdensity(Normal(;μ,τ), y)
        @test ℓ ≈ logdensity(Normal(;μ,logσ), y)
    end
end

@testset "Kernel" begin
    κ = MeasureTheory.kernel(identity, MeasureTheory.Dirac)
    @test rand(κ(1.1)) == 1.1
end

@testset "SpikeMixture" begin
    @test rand(SpikeMixture(Dirac(0), 0.5)) == 0
    @test rand(SpikeMixture(Dirac(1), 1.0)) == 1
    w = 1/3
    m = SpikeMixture(Normal(), w)
    bm = basemeasure(m)
    @test (bm.s*bm.w)*bm.m == 1.0*basemeasure(Normal())
    @test density(m, 1.0)*(bm.s*bm.w) == w*density(Normal(),1.0)
    @test density(m, 0)*(bm.s*(1-bm.w)) ≈ (1-w)
end

@testset "Dirac" begin
    @test rand(Dirac(0.2)) == 0.2
    @test logdensity(Dirac(0.3), 0.3) == 0.0
    @test logdensity(Dirac(0.3), 0.4) == -Inf
end

import MeasureTheory.:⋅
function ⋅(μ::Normal, kernel) 
    m = kernel(μ)
    Normal(μ = m.μ.μ, σ = sqrt(m.μ.σ^2 + m.σ^2))
end

"""
    ConstantMap(β)
Represents a function `f = ConstantMap(β)`
such that `f(x) == β`.
"""
struct ConstantMap{T}
    x::T
end
(a::ConstantMap)(x) = a.x
(a::ConstantMap)() = a.x

struct AffineMap{S,T}
    B::S
    β::T
end
(a::AffineMap)(x) = a.B*x + a.β
(a::AffineMap)(p::Normal) = Normal(μ = a.B*mean(p) + a.β, σ = sqrt(a.B*p.σ^2*a.B'))

@testset "DynamicFor" begin
    mc = Chain(Normal(μ=0.0)) do x Normal(μ=x) end
    @test_broken rand(mc)
    r = rand(Random.GLOBAL_RNG, mc)
   
    # Check that `r` is now deterministic
    @test logdensity(mc, take(r, 100)) == logdensity(mc, take(r, 100))
    
    d2 = For(r) do x Normal(μ=x) end  

    @test_broken let r2 = rand(d2)
        logdensity(d2, take(r2, 100)) == logdensity(d2, take(r2, 100))
    end
end

@testset "Univariate chain" begin
    ξ0 = 1.
    x = 1.2
    P0 = 1.0

    Φ = 0.8
    β = 0.1
    Q = 0.2

    μ = Normal(μ=ξ0, σ=sqrt(P0))
    kernel = MeasureTheory.kernel(Normal; μ=AffineMap(Φ, β), σ=ConstantMap(Q))
    
    @test (μ ⋅ kernel).μ == Normal(μ = 0.9, σ = 0.824621).μ
    
    chain = Chain(kernel, μ)
    

    dyniterate(iter::TimeLift, ::Nothing) = dyniterate(iter, 0=>nothing) 
    tr1 = trace(TimeLift(chain), nothing, u -> u[1] > 15)
    tr2 = trace(TimeLift(rand(Random.GLOBAL_RNG, chain)), nothing, u -> u[1] > 15)
    collect(Iterators.take(chain, 10))
    collect(Iterators.take(rand(Random.GLOBAL_RNG, chain), 10))
    
end
