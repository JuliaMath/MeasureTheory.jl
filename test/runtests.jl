using MeasureTheory
using Test
using StatsFuns
using TransformVariables: transform

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
