using MeasureTheory
using Test

@testset "MeasureTheory.jl" begin
    # Write your own tests here.
end

@testset "Kernel" begin
    κ = MeasureTheory.kernel(MeasureTheory.Dirac, identity)
    @test rand(κ(1.1)) == 1.1
end
