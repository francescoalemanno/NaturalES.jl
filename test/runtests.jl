using NaturalES
using Test
using LinearAlgebra
f(x) = sum(cosh.(x.-(1:length(x)))) - length(x)
g(P) = sum(0<rand()<p<=1 ? 1 : 5+p^2 for p in P)

@testset "Simple Functions" begin
    @test separable_nes(f,(0.0,0.0,0.0)).cost < 1e-5
    @test separable_nes(g,(5.0,5.0,5.0),(1.0,1.0,1.0),0.1,0.01).cost ==3
    @test xnes(f,[0.0,0.0]).cost < 1e-5
end
using BenchmarkTools

xi=zeros(10)
@btime xnes(f,$xi) # 82.308 ms (187816 allocations: 46.74 MiB)
