using NaturalES
using Test

f(x) = sum(cosh.(x.-(1:length(x)))) - length(x)
g(P) = sum(0<rand()<p<=1 ? 1 : 5+p^2 for p in P)

@testset "Simple Functions" begin
    @test separable_nes(f,(0.0,0.0,0.0)).cost < 1e-5
    @test separable_nes(g,(5.0,5.0,5.0),(1.0,1.0,1.0),0.1,0.01).cost ==3
end
#=
f(x)=+(x[1]-2)^2+(x[2]-3)^2
f(x) = sum(cosh.(x.-(1:length(x)))) - length(x)

x=[1.0,2.0]

A=0.0.+1*I(2)|>collect

xnes(f,x,A)
=#
