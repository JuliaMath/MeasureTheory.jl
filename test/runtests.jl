using MeasureTheory
using Test
using StatsFuns

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
        (n,p) = (10,0.3)
        logitp = logit(p)
        probitp = norminvcdf(p)
        (x,y) = draw2(Binomial(n,p))

        let Δlogπ(n,p,x) = logdensity(Binomial(n,p), x) - binomlogpdf(n,p,x)
            @test Δlogπ(n,p,x) ≈ Δlogπ(n,p,y) atol=1e-8
        end

        @test logdensity(Binomial(;n,p), x) ≈ logdensity(Binomial(; n, logitp), x)

        @test logdensity(Binomial(;n,p), x) ≈ logdensity(Binomial(; n, probitp), x)

        @test_broken logdensity(Binomial(n,p), CountingMeasure(ℤ[0:n]), x) ≈ binomlogpdf(n,p,x)
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
