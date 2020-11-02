using MeasureTheory
using Test

@testset "MeasureTheory.jl" begin
    # Write your own tests here.
end

@testset "Kernel" begin
    κ = MeasureTheory.kernel(identity, MeasureTheory.Dirac)
    @test rand(κ(1.1)) == 1.1
end
