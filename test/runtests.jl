using MeasureTheory
using Test

@testset "MeasureTheory.jl" begin
    # Write your own tests here.
end

@testset "Kernel" begin
    κ = MeasureTheory.kernel(x->(x,), MeasureTheory.Dirac)
    @test rand(κ(1.1)) == 1.1
end
