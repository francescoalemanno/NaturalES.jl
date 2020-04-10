module NaturalES
using Random
function optimize(f,x0,σ,lr)
    N=length(x0)
    ϵ=copy(x0)
    ∇f=copy(x0)
    x=copy(x0)
    Z=zero(f(x0))
    samples=4N
    #bestavgF=
    for steps in 1:1000
        avgF=Z
        for i in 1:samples
            randn!(ϵ)
            F=f(x .+ σ.*ϵ)
            avgF+=F
            ∇f.+=F.*ϵ./σ
        end
        avgF/=samples
        ∇f./=samples
        x.-=lr.*∇f
        ∇f.=0
        @show avgF
    end
    x
end
end # module
