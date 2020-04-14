using NaturalES
using Test
using LinearAlgebra

f(x) = sum(cosh.(x.-(1:length(x)))) - length(x)

g(P) = sum(0<rand()<p<=1 ? 1 : 5+p^2 for p in P)

function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

@testset "Simple Functions" begin
    @test separable_nes(f,(0.0,0.0,0.0)).cost < 1e-5
    @test separable_nes(g,(5.0,5.0,5.0),(1.0,1.0,1.0),0.1,0.01).cost ==3
    @test xnes(f,[0.0,0.0]).cost < 1e-5
    @test xnes(rosenbrock2d,[0.0,0.0]).cost < 1e-5
end
