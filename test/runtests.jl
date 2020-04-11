using NaturalES
using Test

f(x) = sum(cosh.(x)) - length(x)

@testset "Simple Functions" begin
    @test optimize(f,[1.0,0.5,2.0],0.1,0.1)[2] < 1e5
    @test optimize(f,[1.0,0.5],0.1,0.1)[2] < 1e5
    @test optimize(f,[1.0 0.5; 0.5 1.0],0.1,0.1)[2] < 1e5
end
