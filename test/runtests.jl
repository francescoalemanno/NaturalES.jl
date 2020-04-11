using NaturalES
using Test

f(x) = sum(cosh.(x.-(1:length(x)))) - length(x)
g(p) = rand()<p[1] ? 1 : 2

@testset "Simple Functions" begin
    @test optimize(f,(1.0,0.5,2.0),0.1,0.1)[2] < 1e-5
    @test optimize(f,(1.0,0.5),0.1,0.1)[2] < 1e-5
    @test optimize(f,(1.0, 0.5, 0.5, 1.0),0.1,0.1)[2] < 1e-5
    @test optimize(g,(0.2,),0.25,0.1)[1][1] > 1
    @test optimize(g,0.2,0.25,0.1)[1] > 1
end
