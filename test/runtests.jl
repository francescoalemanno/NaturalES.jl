using NaturalES
using Test
using LinearAlgebra
using Random
const frng=MersenneTwister(123)

f(x) = sum(cosh.(x.-(1:length(x)))) - length(x)

g(P) = sum(0<rand(frng)<p<=1 ? 1 : 5+p^2 for p in P)

function rosenbrock2d(x)
    s=(1.0 - x[1])^2
    for i in 2:length(x)
        s+=100.0 * (x[i] - x[i-1]^2)^2
    end
    return s
end

@testset "Simple Functions" begin
    @test separable_nes(f,[0.0,0.0,0.0],1.0).cost < 1e-5
    @test separable_nes(f,[4.0,4.0,4.0],1.0).cost < 1e-5
    for N in 2:5
        @test separable_nes(rosenbrock2d,zeros(N),1.0,(ημ=1.0,ησ=0.0005)).cost < 1e-5
    end
    for i in 1:10
        s=separable_nes(g,[2.0,2.0,2.0],1.0,(ημ=0.1,ησ=0.05)).sol
        @test all(0.75 .< s .< 1)
    end
    @test exponential_nes(f,[0.0,0.0],1.0).cost < 1e-5
    @test exponential_nes(rosenbrock2d,[0.5,0.5],1.0).cost < 1e-5
end
