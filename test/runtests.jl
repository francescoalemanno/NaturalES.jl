using NaturalES
using Test
using LinearAlgebra

f(x) = sum(cosh.(x.-(1:length(x)))) - length(x)
sphere(x) = sum(abs2,x)
function rosenbrock2d(x::AbstractVector{T}) where T
    s=(1.0 - x[1])^2
    for i in 1:(length(x)-1)
        s+=100.0 * (x[i+1] - x[i]^2)^2
    end
    return s
end

@testset "Optimizers" begin
    for method in [sNES,xNES], x in [0.5,1.0,4.0], fun in [f,sphere,rosenbrock2d]
        ix=[x,-x]
        @info "method: $method, testing init condition $ix, function: $fun"
        @test optimize(fun,ix,1.0,method).cost < 1e-7
        @test optimize(fun,ix,1.0,method,samples=10).cost < 1e-7
        @test optimize(fun,ix,2.0,method,ησ=0.01,atol=1e-8,iterations=10^4).cost < 1e-7
        @test optimize(fun,ix,2.0,method,ησ=0.01,atol=1e-8,iterations=10).cost > 1e-7
    end
end
