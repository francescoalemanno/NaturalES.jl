module NaturalES
trandn(v::NTuple{N,T}) where {N,T} = ntuple(i->randn(T),Val(N))

function optimize(f,x0::NTuple{N,T},σ::T,lr::T) where {T,N}
    x=x0
    Zv=zero.(x0)
    samples=3N
    MAXWASTEEVAL=500
    improv=MAXWASTEEVAL
    bestavgF=10f(x0)
    Z=zero(bestavgF)
    iσ=inv.(σ)
    while improv > 0
        avgF=Z
        ∇f=Zv
        for i in 1:samples
            ϵ=trandn(x0)
            Fa = f(x .+ σ .* ϵ)
            Fb = f(x .- σ .* ϵ)
            avgF += (Fa + Fb) / 2
            δ = (Fa - Fb) / 2
            ∇f = ∇f .+ δ .* ϵ .* iσ
        end
        avgF /= samples
        if avgF < bestavgF
            improv = MAXWASTEEVAL
            bestavgF = avgF
        else
            improv -= 1
        end
        ∇f = ∇f ./ samples
        x = x .- lr .* ∇f
    end
    (sol=x,cost=f(x))
end
function optimize(f,x0::T,σ::T,lr::T) where T
    @inline g(x)=f(x[1])
    r=optimize(g,(x0,),σ,lr)
    (sol=r.sol[1],cost=r.cost)
end
export optimize
end # module
