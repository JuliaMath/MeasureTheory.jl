using MeasureTheory
using Test

function draw2(μ)
    x = rand(μ)
    y = rand(μ)
    while x == y
        y = rand(μ)
    end
    return (x,y)
end

using Test
using StatsFuns

@testset "Parameterized Measures" begin

    @testset "Binomial" begin
        (n,p) = (10,0.3)
        logit_p = logit(p)
        (x,y) = draw2(Binomial(n,p))

        let Δlogπ(n,p,x) = logdensity(Binomial(n,p), x) - binomlogpdf(n,p,x)
            @test Δlogπ(n,p,x) ≈ Δlogπ(n,p,y)
        end

        @test logdensity(Binomial(;n,p), x) ≈ logdensity(Binomial(; n, logit_p), x)

        @test_broken logdensity(Binomial(n,p), CountingMeasure(ℤ[0:n]), x) ≈ binomlogpdf(n,p,x)
    end

end

@testset "Kernel" begin
    κ = MeasureTheory.kernel(identity, MeasureTheory.Dirac)
    @test rand(κ(1.1)) == 1.1
end
