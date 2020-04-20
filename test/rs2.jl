using NaturalES

sols=[]
cnt=0
function rosenbrock(x::AbstractVector{T}) where T
    global cnt=(cnt+1)%20
    cnt==10 && push!(sols,copy(x))
    s=(1.0 - x[1])^2
    for i in 1:(length(x)-1)
        s+=100.0 * (x[i+1] - x[i]^2)^2
    end
    return s
end

exponential_nes(rosenbrock,[-0.41,0.17],0.01)
sols
A=hcat(sols...)'|>collect

using PyPlot

pygui(true)
function rosenbrock2d(x,y)
  return (1.0 - x)^2 + 100.0 * (y - x^2)^2
end

x=-0.5:0.001:1.5
y=-0.5:0.001:1.5
M=rosenbrock2d.(x,y')
plt.style.use("dark_background")
plot(A[:,2],A[:,1],color="red")
contour(x,y,M,log.(1.01:0.5:10))
colorbar()
scatter(1,1,color="blue")


using NaturalES
function rosenbrock2d(x::AbstractVector{T}) where T
    s=(1.0 - x[1])^2
    for i in 1:(length(x)-1)
        s+=100.0 * (x[i+1] - x[i]^2)^2
    end
    return s
end

using BenchmarkTools
@btime exponential_nes($rosenbrock2d,$[0.0,0.0],$1.0)
