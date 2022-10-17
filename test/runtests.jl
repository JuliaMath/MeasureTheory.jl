using Test
using StatsFuns
using Base.Iterators: take
using Random
using LinearAlgebra
# using DynamicIterators: trace, TimeLift
using TransformVariables: transform, as𝕀
using FillArrays

using MeasureTheory
using MeasureBase.Interface
using MeasureTheory: kernel
using Aqua
using IfElse

# Aqua._test_ambiguities(
#     Aqua.aspkgids(MeasureTheory);
#     exclude = [Random.AbstractRNG],
#     # packages::Vector{PkgId};
#     # color::Union{Bool, Nothing} = nothing,
#     # exclude::AbstractArray = [],
#     # # Options to be passed to `Test.detect_ambiguities`:
#     # detect_ambiguities_options...,
# )

Aqua.test_all(MeasureBase; ambiguities = false)

function draw2(μ)
    x = rand(μ)
    y = rand(μ)
    while x == y
        y = rand(μ)
    end
    return (x, y)
end

x = randn(10, 3)
Σ = cholesky(x' * x)
Λ = cholesky(inv(Σ))
σ = MeasureTheory.getL(Σ)
λ = MeasureTheory.getL(Λ)

test_measures = Any[
    # Chain(x -> Normal(μ=x), Normal(μ=0.0))
    For(3) do j
        Normal(σ = j)
    end
    For(2, 3) do i, j
        Normal(i, j)
    end
    Normal()^3
    Normal()^(2, 3)
    3 * Normal()
    Bernoulli(0.2)
    Beta(2, 3)
    Binomial(10, 0.3)
    BetaBinomial(10, 2, 3)
    Cauchy()
    Dirichlet(ones(3))
    Exponential()
    Gumbel()
    Laplace()
    LKJCholesky(3, 2.0)
    Multinomial(n = 10, p = [0.2, 0.3, 0.5])
    NegativeBinomial(5, 0.2)
    Normal(2, 3)
    Poisson(3.1)
    StudentT(ν = 2.1)
    MvNormal(σ = [1 0; 0 1; 1 1])
    MvNormal(λ = [1 0 1; 0 1 1])
    MvNormal(Σ = Σ)
    MvNormal(Λ = Λ)
    MvNormal(σ = σ)
    MvNormal(λ = λ)
    Uniform()
    Counting(Float64)
    Dirac(0.0) + Normal()
]

testbroken_measures = Any[
    Pushforward(as𝕀, Normal())
    # InverseGamma(2) # Not defined yet
    # MvNormal(I(3)) # Entirely broken for now
    TrivialMeasure()
]

@testset "testvalue" begin
    for μ in test_measures
        @info "testing $μ"
        test_interface(μ)
    end

    for μ in testbroken_measures
        @info "testing $μ"
        @test_broken test_measure(μ)
    end

    @testset "testvalue(::Chain)" begin
        mc = Chain(x -> Normal(μ = x), Normal(μ = 0.0))
        r = testvalue(mc)
        @test logdensity_def(mc, Iterators.take(r, 10)) isa AbstractFloat
    end
end

@testset "Parameterized Measures" begin
    @testset "Binomial" begin
        D = Binomial{(:n, :p)}
        par = merge((n = 20,), transform(asparams(D, (n = 20,)), randn(1)))
        d = D(par)
        (n, p) = (par.n, par.p)
        logitp = logit(p)
        probitp = norminvcdf(p)
        y = rand(d)

        ℓ = logdensity_def(Binomial(; n, p), y)
        @test ℓ ≈ logdensity_def(Binomial(; n, logitp), y)
        @test ℓ ≈ logdensity_def(Binomial(; n, probitp), y)

        @test_broken logdensity_def(Binomial(n, p), CountingBase(ℤ[0:n]), x) ≈
                     binomlogpdf(n, p, x)
    end

    @testset "Exponential" begin
        r = rand(MersenneTwister(123), Exponential(2))
        @test r ≈ rand(MersenneTwister(123), Exponential(β = 2))
        @test r ≈ rand(MersenneTwister(123), Exponential(λ = 0.5))
        @test r ≈ rand(MersenneTwister(123), Exponential(logβ = log(2)))
        @test r ≈ rand(MersenneTwister(123), Exponential(logλ = log(0.5)))

        ℓ = logdensity_def(Exponential(2), r)
        @test ℓ ≈ logdensity_def(Exponential(β = 2), r)
        @test ℓ ≈ logdensity_def(Exponential(λ = 0.5), r)
        @test ℓ ≈ logdensity_def(Exponential(logβ = log(2)), r)
        @test ℓ ≈ logdensity_def(Exponential(logλ = log(0.5)), r)
    end

    @testset "NegativeBinomial" begin
        D = NegativeBinomial{(:r, :p)}
        par = transform(asparams(D), randn(2))
        d = D(par)
        (r, p) = (par.r, par.p)
        logitp = logit(p)
        λ = p * r / (1 - p)
        logλ = log(λ)
        y = rand(d)

        ℓ = logdensity_def(NegativeBinomial(; r, p), y)
        @test ℓ ≈ logdensity_def(NegativeBinomial(; r, logitp), y)
        @test ℓ ≈ logdensity_def(NegativeBinomial(; r, λ), y)
        @test ℓ ≈ logdensity_def(NegativeBinomial(; r, logλ), y)

        sample1 = rand(MersenneTwister(123), NegativeBinomial(; r, λ))
        sample2 = rand(MersenneTwister(123), NegativeBinomial(; r, logλ))
        @test sample1 == sample2

        @test_broken logdensity_def(Binomial(n, p), CountingBase(ℤ[0:n]), x) ≈
                     binomlogpdf(n, p, x)
    end

    @testset "Poisson" begin
        sample1 = rand(MersenneTwister(123), Poisson(; logλ = log(100)))
        sample2 = rand(MersenneTwister(123), Poisson(; λ = 100))
        @test sample1 == sample2
    end

    # Fails because we need `asparams` for `::AffinePushfwd`
    # @testset "Normal" begin
    #     D = affine{(:μ,:σ), Normal}
    #     par = transform(asparams(D), randn(2))
    #     d = D(par)
    #     @test params(d) == par

    #     μ = par.μ
    #     σ = par.σ
    #     σ² = σ^2
    #     τ = 1/σ²
    #     logσ = log(σ)
    #     y = rand(d)

    #     ℓ = logdensity_def(Normal(;μ,σ), y)
    #     @test ℓ ≈ logdensity_def(Normal(;μ,σ²), y)
    #     @test ℓ ≈ logdensity_def(Normal(;μ,τ), y)
    #     @test ℓ ≈ logdensity_def(Normal(;μ,logσ), y)
    # end

    @testset "LKJCholesky" begin
        D = LKJCholesky{(:k, :η)}
        par = transform(asparams(D, (k = 4,)), randn(1))
        d = D(merge((k = 4,), par))
        # @test params(d) == par

        η = par.η
        logη = log(η)

        y = rand(d)
        η = par.η
        ℓ = logdensity_def(LKJCholesky(4, η), y)
        @test ℓ ≈ logdensity_def(LKJCholesky(k = 4, logη = logη), y)
    end
end

@testset "TransitionKernel" begin
    κ = kernel(Dirac)
    @test rand(κ(1.1)) == 1.1

    k1 = kernel() do x
        Normal(x, x^2)
    end

    k2 = kernel(Normal) do x
        (μ = x, σ = x^2)
    end

    k3 = kernel(Normal; μ = identity, σ = abs2)

    k4 = kernel(Normal; μ = first, σ = last) do x
        (x, x^2)
    end

    @test k1(3) == k2(3) == k3(3) == k4(3) == Normal(3, 9)
end

@testset "For" begin
    FORDISTS = [
        For(1:10) do j
            Normal(μ = j)
        end
        For(4, 3) do μ, σ
            Normal(μ, σ)
        end
        For(1:4, 1:4) do μ, σ
            Normal(μ, σ)
        end
        For(eachrow(rand(4, 2))) do x
            Normal(x[1], x[2])
        end
        For(rand(4), rand(4)) do μ, σ
            Normal(μ, σ)
        end
    ]

    for d in FORDISTS
        @info "testing $d"
        @test logdensity_def(d, rand(d)) isa Float64
    end
end

import MeasureTheory: ⋅

function ⋅(μ::Normal, kernel)
    m = kernel(μ)
    Normal(μ = m.μ.μ, σ = sqrt(m.μ.σ^2 + m.σ^2))
end
struct AffinePushfwdMap{S,T}
    B::S
    β::T
end
(a::AffinePushfwdMap)(x) = a.B * x + a.β
function (a::AffinePushfwdMap)(p::Normal)
    Normal(μ = a.B * mean(p) + a.β, σ = sqrt(a.B * p.σ^2 * a.B'))
end

@testset "DynamicFor" begin
    mc = Chain(Normal(μ = 0.0)) do x
        Normal(μ = x)
    end
    r = rand(mc)

    # Check that `r` is now deterministic
    @test logdensity_def(mc, take(r, 100)) == logdensity_def(mc, take(r, 100))

    d2 = For(r) do x
        Normal(μ = x)
    end

    @test let r2 = rand(d2)
        logdensity_def(d2, take(r2, 100)) == logdensity_def(d2, take(r2, 100))
    end
end

@testset "Product of Diracs" begin
    x = randn(3)
    t = as(productmeasure(Dirac.(x)))
    @test transform(t, []) == x
end

# @testset "Univariate chain" begin
#     ξ0 = 1.
#     x = 1.2
#     P0 = 1.0

#     Φ = 0.8
#     β = 0.1
#     Q = 0.2

#     μ = Normal(μ=ξ0, σ=sqrt(P0))
#     kernel = MeasureTheory.kernel(Normal; μ=AffinePushfwdMap(Φ, β), σ=MeasureTheory.AsConst(Q))

#     @test (μ ⋅ kernel).μ == Normal(μ = 0.9, σ = 0.824621).μ

#     chain = Chain(kernel, μ)

#     dyniterate(iter::TimeLift, ::Nothing) = dyniterate(iter, 0=>nothing)
#     tr1 = trace(TimeLift(chain), nothing, u -> u[1] > 15)
#     tr2 = trace(TimeLift(rand(Random.GLOBAL_RNG, chain)), nothing, u -> u[1] > 15)
#     collect(Iterators.take(chain, 10))
#     collect(Iterators.take(rand(Random.GLOBAL_RNG, chain), 10))
# end

@testset "rootmeasure/logpdf" begin
    x = rand(Normal())
    @test logdensityof(𝒹(Normal(), rootmeasure(Normal())), x) ≈ logdensityof(Normal(), x)
end

@testset "Transforms" begin
    t = as𝕀
    @testset "Pushforward" begin
        μ = Normal()
        ν = Pushforward(t, μ)
        x = rand(μ)
        @test logdensity_def(μ, x) ≈ logdensity_def(Pushforward(TV.inverse(t), ν), x)
    end

    @testset "Pullback" begin
        ν = Uniform()
        μ = Pullback(t, ν)
        y = rand(ν)
        @test logdensity_def(ν, y) ≈ logdensity_def(Pullback(TV.inverse(t), μ), y)
    end
end

@testset "Likelihood" begin
    dps = [(Normal(), 2.0)
    # (Pushforward(as((μ=asℝ,)), Normal()^1), (μ=2.0,))
    ]

    ℓs = [
        Likelihood(Normal{(:μ,)}, 3.0),
        # Likelihood(kernel(Normal, x -> (μ=x, σ=2.0)), 3.0)
    ]

    for (d, p) in dps
        for ℓ in ℓs
            @test logdensityof(d ⊙ ℓ, p) ≈ logdensityof(d, p) + logdensityof(ℓ.k(p), ℓ.x)
        end
    end
end

@testset "Reproducibility" begin

    # NOTE: The `test_broken` below are mostly because of the change to `AffinePushfwd`.
    # For example, `Normal{(:μ,:σ)}` is now `AffinePushfwd{(:μ,:σ), Normal{()}}`.
    # The problem is not really with these measures, but with the tests
    # themselves.
    #
    # We should instead probably be doing e.g.
    # `D = typeof(Normal(μ=0.3, σ=4.1))`

    function repro(D, args, nt = NamedTuple())
        t = asparams(D{args}, nt)
        d = D(transform(t, randn(t.dimension)))
        r(d) = rand(Random.MersenneTwister(1), d)
        logdensity_def(d, r(d)) == logdensity_def(d, r(d))
    end

    @testset "Bernoulli" begin
        @test repro(Bernoulli, (:p,))
    end
    @testset "Binomial" begin
        @test repro(Binomial, (:n, :p), (n = 10,))
    end

    @testset "Beta" begin
        @test repro(Beta, (:α, :β))
    end

    @testset "BetaBinomial" begin
        @test repro(BetaBinomial, (:n, :α, :β), (n = 10,))
    end

    @testset "Cauchy" begin
        @test_broken repro(Cauchy, (:μ, :σ))
    end

    @testset "Dirichlet" begin
        @test_broken repro(Dirichlet, (:p,))
    end

    @testset "Exponential" begin
        @test repro(Exponential, (:λ,))
    end

    @testset "Gumbel" begin
        @test_broken repro(Gumbel, (:μ, :σ))
    end

    @testset "InverseGamma" begin
        @test_broken repro(InverseGamma, (:p,))
    end

    @testset "Laplace" begin
        @test_broken repro(Laplace, (:μ, :σ))
    end

    @testset "LKJCholesky" begin
        @test repro(LKJCholesky, (:k, :η), (k = 3,))
    end

    @testset "Multinomial" begin
        @test_broken repro(Multinomial, (:n, :p))
    end

    @testset "MvNormal" begin
        @test_broken repro(MvNormal, (:μ,))

        σ = LowerTriangular(randn(3, 3))
        Σ = σ * σ'
        d = MvNormal(σ = σ)
        x = rand(d)
        @test logdensityof(d, x) ≈ logdensityof(Dists.MvNormal(Σ), x)
        @test logdensityof(MvNormal(zeros(3), σ), x) ≈ logdensityof(d, x)
    end

    @testset "NegativeBinomial" begin
        @test repro(NegativeBinomial, (:r, :p))
    end

    @testset "Normal" begin
        @test_broken repro(Normal, (:μ, :σ))
    end

    @testset "Poisson" begin
        @test repro(Poisson, (:λ,))
    end

    @testset "StudentT" begin
        @test_broken repro(StudentT, (:ν, :μ))
    end

    @testset "Uniform" begin
        @test repro(Uniform, ())
    end
end

@testset "ProductMeasure" begin
    d = For(1:10) do j
        Poisson(exp(j))
    end
    x = Vector{Int16}(undef, 10)
    @test rand!(d, x) isa Vector
    @test rand(d) isa Vector

    @testset "Indexed by Generator" begin
        d = For((j^2 for j in 1:10)) do i
            Poisson(i)
        end
        x = Vector{Int16}(undef, 10)
        @test rand!(d, x) isa Vector
        @test rand(d) isa MeasureTheory.RealizedSamples
    end

    @testset "Indexed by multiple Ints" begin
        d = For(2, 3) do μ, σ
            Normal(μ, σ)
        end
        x = Matrix{Float16}(undef, 2, 3)
        @test rand!(d, x) isa Matrix
        @test_broken rand(d) isa Matrix{Float16}
    end
end

@testset "Show methods" begin
    @testset "PowerMeasure" begin
        @test repr(Lebesgue(ℝ)^5) == "Lebesgue(ℝ) ^ 5"
        @test repr(Lebesgue(ℝ)^(3, 2)) == "Lebesgue(ℝ) ^ (3, 2)"
    end
end

@testset "Density measures and Radon-Nikodym" begin
    x = randn()
    let d = ∫(𝒹(Cauchy(), Normal()), Normal())
        @test logdensityof(𝒹(d, Cauchy()), x) ≈ 0 atol = 1e-12
    end

    let f = 𝒹(∫(x -> x^2, Normal()), Normal())
        @test densityof(f, x) ≈ x^2
    end

    # let d = ∫exp(log𝒹(Cauchy(), Normal()), Normal())
    #     @test logdensity_def(d, Cauchy(), x) ≈ 0 atol=1e-12
    # end

    let f = 𝒹(∫exp(x -> x^2, Normal()), Normal())
        @test logdensityof(f, x) ≈ x^2
    end
end

@testset "Half measures" begin
    @testset "HalfNormal" begin
        d = Normal(σ = 3)
        h = HalfNormal(3)
        x = rand(h)
        @test densityof(𝒹(h, Lebesgue(ℝ)), x) ≈ 2 * densityof(𝒹(d, Lebesgue(ℝ)), x)
    end

    @testset "HalfCauchy" begin
        d = Cauchy(σ = 3)
        h = HalfCauchy(3)
        x = rand(h)
        @test densityof(𝒹(h, Lebesgue(ℝ)), x) ≈ 2 * densityof(𝒹(d, Lebesgue(ℝ)), x)
    end

    @testset "HalfStudentT" begin
        d = StudentT(ν = 2, σ = 3)
        h = HalfStudentT(2, 3)
        x = rand(h)
        @test densityof(𝒹(h, Lebesgue(ℝ)), x) ≈ 2 * densityof(𝒹(d, Lebesgue(ℝ)), x)
    end
end

@testset "MvNormal" begin
    Q, R = qr(randn(4, 2))
    D = Diagonal(sign.(diag(R)))
    Q = Matrix(Q) * D
    R = D * R

    z = randn(2)

    ℓ = logdensityof(MvNormal((σ = R,)), z)
    @test ℓ ≈ logdensityof(Dists.MvNormal(R * R'), z)
    @test ℓ ≈ logdensityof(MvNormal((σ = Q * R,)), Q * z)
end

@testset "AffinePushfwd" begin
    testmeasures = [
        (Normal, NamedTuple())
        (Cauchy, NamedTuple())
        (Laplace, NamedTuple())
        (StudentT, (ν = 3,))
    ]

    using Test
    function test_noerrors(d)
        x = rand(d)
        @test logdensityof(d, x) isa Real
        @test logdensity_def(d, x) isa Real
    end

    for (M, nt) in testmeasures
        for p in [(μ = 1,), (μ = 1, σ = 2), (μ = 1, λ = 2), (σ = 2,), (λ = 2,)]
            d = M(merge(nt, p))
            @info "Testing $d"
            test_noerrors(d)
        end
        # for n in 1:3
        #     @show n
        #     for k in 1:n
        #         @show k
        #         pars = [(μ=randn(n),), (μ=randn(n),σ=randn(n,k)), (μ=randn(n),λ=randn(k,n)), (σ=randn(n,k),), (λ=randn(k,n),)]
        #         for p in pars
        #             @show p
        #             d = M(merge(nt, p))
        #             test_noerrors(d)
        #         end
        #     end
        # end
    end
end

@testset "AffineTransform" begin
    f = AffineTransform((μ = 3, σ = 2))
    @test f(inverse(f)(1)) == 1
    @test inverse(f)(f(1)) == 1

    f = AffineTransform((μ = 3, λ = 2))
    @test f(inverse(f)(1)) == 1
    @test inverse(f)(f(1)) == 1

    f = AffineTransform((σ = 2,))
    @test f(inverse(f)(1)) == 1
    @test inverse(f)(f(1)) == 1

    f = AffineTransform((λ = 2,))
    @test f(inverse(f)(1)) == 1
    @test inverse(f)(f(1)) == 1

    f = AffineTransform((μ = 3,))
    @test f(inverse(f)(1)) == 1
    @test inverse(f)(f(1)) == 1

    f = AffineTransform((σ = [1 2; 2 1],))
    @test f(inverse(f)([1 2; 2 1])) == [1 2; 2 1]
    @test inverse(f)(f([1 2; 2 1])) == [1 2; 2 1]
end

@testset "AffinePushfwd" begin
    unif = ∫(x -> 0 < x < 1, Lebesgue(ℝ))
    f1 = AffineTransform((μ = 3.0, σ = 2.0))
    f2 = AffineTransform((μ = 3.0, λ = 2.0))
    f3 = AffineTransform((μ = 3.0,))
    f4 = AffineTransform((σ = 2.0,))
    f5 = AffineTransform((λ = 2.0,))

    for f in [f1, f2, f3, f4, f5]
        par = getfield(f, :par)
        @test AffinePushfwd(par)(unif) == AffinePushfwd(f, unif)
        @test densityof(AffinePushfwd(f, AffinePushfwd(inverse(f), unif)), 0.5) == 1
    end

    d = ∫exp(x -> -x^2, Lebesgue(ℝ))

    μ = randn(3)
    # σ = LowerTriangular(randn(3, 3))
    σ = let x = randn(10, 3)
        cholesky(x' * x).L
    end
    λ = inv(σ)

    x = randn(3)

    @test logdensity_def(AffinePushfwd((μ = μ, σ = σ), d^3), x) ≈
          logdensity_def(AffinePushfwd((μ = μ, λ = λ), d^3), x)
    @test logdensity_def(AffinePushfwd((σ = σ,), d^3), x) ≈
          logdensity_def(AffinePushfwd((λ = λ,), d^3), x)
    @test logdensity_def(AffinePushfwd((μ = μ,), d^3), x) ≈ logdensity_def(d^3, x - μ)

    d = ∫exp(x -> -x^2, Lebesgue(ℝ))
    a = AffinePushfwd((σ = [1 0]',), d^1)
    x = randn(2)
    y = randn(1)
    @test logdensityof(a, x) ≈ logdensityof(d, inverse(a.f)(x)[1])
    @test logdensityof(a, a.f(y)) ≈ logdensityof(d^1, y)

    b = AffinePushfwd((λ = [1 0]'',), d^1)
    @test logdensityof(b, x) ≈ logdensityof(d, inverse(b.f)(x)[1])
    @test logdensityof(b, b.f(y)) ≈ logdensityof(d^1, y)
end

@testset "IfElseMeasure" begin
    p = rand()
    x = randn()

    @test let
        a = logdensityof(IfElse.ifelse(Bernoulli(p), Normal(), Normal()), x)
        b = logdensityof(Normal(), x)
        a ≈ b
    end

    @test let
        a = logdensityof(IfElse.ifelse(Bernoulli(p), Normal(2, 3), Normal()), x)
        b = logdensityof(p * Normal(2, 3) + (1 - p) * Normal(), x)
        a ≈ b
    end
end

@testset "https://github.com/cscherrer/MeasureTheory.jl/issues/217" begin
    d = For(rand(3), rand(3)) do x, y
        Normal(x, y)
    end

    x = rand(d)

    @test logdensityof(d, x) isa Real
end

@testset "Distributions.jl cdf" begin
    @test cdf(Normal(0, 1), 0) == 0.5
    @test cdf.((Normal(0, 1),), [0, 0]) == [0.5, 0.5]
end
