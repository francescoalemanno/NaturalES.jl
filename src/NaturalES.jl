module NaturalES
using Random
using LinearAlgebra
include("symprod.jl")
trandn(rn::AbstractRNG, v::NTuple{N,T}) where {N,T} = ntuple(i->randn(rn,T),Val(N))
function tsqnorm(v::NTuple{N,T}) where {N,T}
    s=zero(T)
    @fastmath for i in eachindex(v)
        e=@inbounds v[i]
        s+=e*e
    end
    s
end

function separable_nes(f,x0::NTuple{N,T},σ::NTuple{N,T},lr_x::T,lr_σ::T) where {T,N}
    rn=MersenneTwister(123+3N)
    x=x0
    Zv=zero.(x0)
    samples=4 + ceil(Int, log(3*N))
    MAXWASTEEVAL=500*samples
    improv=MAXWASTEEVAL
    bestavgF=10f(x0)
    Z=zero(bestavgF)
    O=one(T)
    mirror=(-O,O)
    Nm=length(mirror)
    while improv > 0
        avgF=Z
        ∇f=Zv
        ∇fσ=Zv
        for i in 1:samples
            ϵ=trandn(rn,x0)
            for s in mirror
                F = f(x .+ s .* σ .* ϵ)
                avgF += F
                ∇f = ∇f .+ s*F .* ϵ
                ∇fσ = ∇fσ .+ F .* (ϵ.*ϵ .- O)
            end
        end
        avgF /= Nm*samples
        if avgF < bestavgF
            improv = MAXWASTEEVAL
            bestavgF = avgF
        else
            improv -= 1
        end
        g=(∇fσ...,∇f...)

        if !all(isfinite,g)
            throw(ErrorException("Function produced $avgF near $x :\n try lowering σ, the learning rate, or changing initial point."))
        end

        ∇f = clamp.(∇f ./ (Nm*samples),-O,O)
        ∇fσ = clamp.(∇fσ ./ (Nm*samples),-O/lr_σ,O/lr_σ)
        x = x .- lr_x .* σ .* ∇f
        σ = σ .* exp.(-lr_σ ./ 2 .*∇fσ)
        #@show x,σ,avgF
    end
    (sol=x, cost=f(x))
end

function separable_nes(f,x0::NTuple{N,T},
    σ::NTuple{N,T}=ntuple(i->one(T),Val(N)),
    lr_x=one(T),
    lr_σ=convert(T,(3 + log(N)) / (5 * sqrt(N)))
    ) where {N,T}
    separable_nes(f,x0,vσ,lr_x,lr_σ)
end

function xnes(f,μ::AbstractVector{T},A::AbstractMatrix{T}) where T
    d=length(μ)
    n=4+floor(Int,3*log(d))
    Z=[Array{T}(undef,d) for i in 1:n]
    SZ=SymProd.(Z)
    σ=T(abs(det(A))^(1/d))
    F=fill(f(μ),n)
    B=A./σ
    idx=collect(1:n)
    u=@. max(0,log(n/2+1)-log(idx))
    norm_u=sum(u)
    @. u=T(u/norm_u - 1/n)
    ημ=one(T)
    ησ=ηB=T(3/5*(3+log(d))/(d*√d))
    Gδ=zeros(T,d)
    GM=zeros(T,d,d)
    GB=zeros(T,d,d)
    Id=I(d)
    while σ>1e-5
        for i in 1:n
            randn!(Z[i])
            F[i]=f(μ .+ σ.*B*Z[i])
        end
        sort!(idx,by=i->F[i])
        @. begin
            Gδ=u[1]*Z[idx[1]]
            GM=u[1]*(SZ[idx[1]]-Id)
            for i in 2:n
                Gδ+=u[i]*Z[idx[i]]
                GM+=u[i]*(SZ[idx[i]]-Id)
            end
        end
        Gσ=tr(GM)/d
        GM .-= Gσ
        μ.+=ημ*σ*B*Gδ
        σ=σ*exp(ησ/2*Gσ)
        B.*=exp(ηB/2 .* GM)
    end
    return (sol=μ, cost=f(μ))
end

function xnes(f,μ::AbstractVector{T},σ::AbstractVector{T}) where T
    A=diagm(σ)
    xnes(f,μ,A)
end

function xnes(f,μ::AbstractVector{T},σ::T=one(T)) where T
    xnes(f,μ,fill(σ,length(μ)))
end

export separable_nes,xnes
end # module
