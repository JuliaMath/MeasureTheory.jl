using Test
using StatsFuns
using Base.Iterators: take
using Random
using LinearAlgebra
using DynamicIterators: trace, TimeLift
using TransformVariables: transform, as𝕀, inverse

if Base.VERSION ≥ v"1.6"
    @testset "No warnings on import" begin
        @test_nowarn @eval using MeasureTheory
    end
else
    @eval using MeasureTheory
end

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

    @testset "NegativeBinomial" begin
        D = NegativeBinomial{(:r, :p)}
        par = transform(asparams(D), randn(2))
        d = D(par)
        (r,p) = (par.r, par.p)
        logitp = logit(p)
        λ = p * r / (1 - p)
        y = rand(d)

        ℓ = logdensity(NegativeBinomial(;r, p), y)
        @test ℓ ≈ logdensity(NegativeBinomial(;r, logitp), y)
        @test ℓ ≈ logdensity(NegativeBinomial(;r, λ), y)

        @test_broken logdensity(Binomial(n,p), CountingMeasure(ℤ[0:n]), x) ≈ binomlogpdf(n,p,x)
    end

    @testset "Poisson" begin
        sample1 = rand(MersenneTwister(123), Poisson{(:logλ,)}(log(100)))
        sample2 = rand(MersenneTwister(123), Poisson{(:λ,)}(100))
        @test sample1 == sample2
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

    @testset "LKJCholesky" begin
        D = LKJCholesky{4}{(:η,)}
        par = transform(asparams(D), randn(1))
        d = D(par)
        @test params(d) == par

        η  = par.η
        logη = log(η)

        y = rand(d)
        η = par.η
        ℓ = logdensity(LKJCholesky{4}(η), y)
        @test ℓ ≈ logdensity(LKJCholesky{4}(logη=logη), y)
    end
end

@testset "Kernel" begin
    κ = MeasureTheory.kernel(MeasureTheory.Dirac, identity)
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

@testset "For" begin
    FORDISTS = [
        For(1:10) do j Normal(μ=j) end
        For(4,3) do μ,σ Normal(μ,σ) end
        For(1:4, 1:4) do μ,σ Normal(μ,σ) end
        For(eachrow(rand(4,2))) do x Normal(x[1], x[2]) end
        For(rand(4), rand(4)) do μ,σ Normal(μ,σ) end
    ]

    for d in FORDISTS
        @test logdensity(d, rand(d)) isa Float64
    end
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
    r = rand(mc)
   
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

@testset "Transforms" begin
    t = as𝕀
    @testset "Pushforward" begin
        μ = Normal()
        ν = Pushforward(t, μ)
        x = rand(μ)
        @test logdensity(μ, x) ≈ logdensity(Pushforward(inverse(t), ν), x)
    end

    @testset "Pullback" begin
        ν = Uniform()
        μ = Pullback(t,ν)
        y = rand(ν)
        @test logdensity(ν, y) ≈ logdensity(Pullback(inverse(t), μ), y)
    end
end

using TransformVariables

@testset "Likelihood" begin
    dps = [
        (Normal()                             ,    2.0  )
        # (Pushforward(as((μ=asℝ,)), Normal()^1), (μ=2.0,))
    ]

    ℓs = [
        Likelihood(Normal{(:μ,)},              3.0)
        Likelihood(kernel(Normal, x -> (μ=x, σ=2.0)), 3.0)
    ]

    for (d,p) in dps
        for ℓ in ℓs
            @test logdensity(d ⊙ ℓ, p) == logdensity(d, p) + logdensity(ℓ, p)
        end
    end
end
