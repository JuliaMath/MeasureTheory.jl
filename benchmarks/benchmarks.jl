using MeasureTheory, Static, BenchmarkTools, LinearAlgebra, StaticArrays

x = randn(200,200);
Σ = x' * x;

using PositiveFactorizations
CΣ = cholesky(Positive, Σ)


function benchmark()
    σtimes = Float64[]
    ωtimes = Float64[]
    y = randn(200);
    for k in 1:100
        U = UpperTriangular(SizedMatrix{k,k}(CΣ.U.data[1:k,1:k]));
        C = Cholesky(U, 'U', 0);
        yk = y[1:k]
        push!(σtimes, @belapsed logdensity(MvNormal(σ=$C), $yk))
        push!(ωtimes, @belapsed logdensity(MvNormal(ω=$C), $yk))
    end

    (σtimes, ωtimes)
end

(σtimes, ωtimes) = benchmark()

# using PDMats
# C = cholesky(Σ[1:20,1:20])
# Dists.MvNormal(PDMat(C))

# @btime Dists.logpdf(Dists.MvNormalCanon(PDMat($C)), $y)
# @btime Dists.logpdf(Dists.MvNormal(PDMat($C)), $y)

plot(σtimes[1:end] .* 1e6; dpi=300, legend=150, label="σ parameterization")
plot!(ωtimes[1:end] .* 1e6, label="ω parameterization")
xlabel!("Dimensions")
ylabel!("Time (μs)")
title!("MeasureTheory.jl Log-density for MvNormal")
savefig("mvnormal.svg")
savefig("mvnormal.png")

# C_unsized = Cholesky(C.UL.data.data, 'U', )

# using PDMats
# PD = PDMat(C_unsized)

# julia> @btime Dists.logpdf(Dists.MvNormalCanon($PD), $y)
#   4.944 μs (6 allocations: 1.97 KiB)
# -124.71152557023147

# julia> @btime Dists.logpdf(Dists.MvNormalCanon(PDMat($C_unsized)), $y)
#   16.552 μs (10 allocations: 27.25 KiB)
# -124.71152557023147

# julia> @btime logdensity(MvNormal(ω=$C), $y)
#   214.948 ns (0 allocations: 0 bytes)
# -124.53160601776467


# L = inv(C.L)
# Σ = L' * L



# @btime Dists.logpdf(Dists.MvNormalCanon($PD), $y)
# @btime Dists.logpdf(Dists.MvNormalCanon(PDMat($C_unsized)), $y)



Dists.logpdf(Dists.MvNormal(Σ), y)
logdensity(MvNormal(ω=C), y)



# n = 40;
# t = CorrCholesky(StaticInt(n));
# C = transform(t, randn(dimension(t)));
@btime MeasureTheory.logdet_pos($C)

# @btime logdensity(MvNormal(ω=$C), $(randn(n)))
