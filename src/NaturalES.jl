module NaturalES
using Random
function optimize(f,x0::AbstractArray{T},σ,lr) where T
    N=length(x0)
    ϵ=copy(x0)
    ∇f=copy(x0)
    x=copy(x0)
    Z=zero(T)
    samples=3N
    MAXWASTEEVAL=500
    improv=MAXWASTEEVAL
    bestavgF=typemax(T)
    iσ=inv.(σ)
    while improv > 0
        avgF=Z
        for i in 1:samples
            randn!(ϵ)
            Fa = f(x .+ σ .* ϵ)
            Fb = f(x .- σ .* ϵ)
            avgF += (Fa + Fb) / 2
            δ = (Fa - Fb) / 2
            ∇f .+= δ .* ϵ .* iσ
        end
        avgF /= samples
        if avgF < bestavgF
            improv = MAXWASTEEVAL
            bestavgF = avgF
        else
            improv -= 1
        end
        ∇f ./= samples
        x .-= lr .* ∇f
        ∇f .= 0
    end
    (sol=x,cost=f(x))
end
end # module
