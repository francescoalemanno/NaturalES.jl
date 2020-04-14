module NaturalES
using Random
using LinearAlgebra


function separable_nes(f,x0::AbstractVector{T},σ::AbstractVector{T},lr_x::T,lr_σ::T) where T
    N=length(x0)
    length(σ)==N || error("the length of 'σ' must be equal to the length of 'x'.")
    rn=MersenneTwister(123+3N)
    x=copy(x0)
    samples=4 + ceil(Int, log(3*N))
    Z=zero(T)
    O=one(T)
    mirror=(-O,O)
    Nm=length(mirror)
    F=fill(f(x0),samples)
    ∇f=fill(Z,N)
    ∇fσ=fill(Z,N)
    ϵ=map(i->zeros(T,N),1:samples)
    u=@. max(0,log(samples/2+1)-log(1:samples))
    u=T.(u./sum(u) .- 1/samples)
    idx=collect(1:samples)
    tmp_x=copy(x0)
    for i in 1:samples
        randn!(rn,ϵ[i])
        tmp_x .= x .+ σ .* ϵ[i]
        F[i] = f(tmp_x)
    end
    sort!(idx,by=i->F[i])
    while sum(σ) > 1e-5
        for i in 1:samples
            j=idx[i]
            if i>=(samples÷2)
                randn!(rn,ϵ[j])
            end
            tmp_x .= x .+ σ .* ϵ[j]
            F[j] = f(tmp_x)
        end
        sort!(idx,by=i->F[i])
        ∇f.=Z
        ∇fσ.=Z
        for i in 1:samples
            j=idx[i]
            ∇f .+= u[i] .* ϵ[j]
            ∇fσ .+= u[i] .* (ϵ[j] .* ϵ[j] .- O)
        end
        x .+= lr_x .* σ .* ∇f
        σ .*= exp.(lr_σ ./ 2 .*∇fσ)
        #@show x,σ,avgF
    end
    (sol=x, cost=f(x))
end

function xnes(f,μ::AbstractVector{T},A::AbstractMatrix{T}) where T
    d=length(μ)
    n=4+floor(Int,3*log(d))
    Z=Array{T}(undef,d,n)
    σ=abs(det(A))^(1/d)
    F=fill(f(μ),n)
    B=A./σ
    idx=collect(1:n)
    u=@. max(0,log(n/2+1)-log(idx))
    u=u./sum(u) .- 1/n
    ημ=one(T)
    ησ=ηB=3/5*(3+log(d))/(d*√d)
    while σ>1e-5
        randn!(Z)
        for i in 1:n
            F[i]=f(μ .+ σ.*B*Z[:,i])
        end
        su=u[invperm(sort(idx,by=i->F[i]))]
        Gδ=sum(su[i].*Z[:,i] for i in 1:n)
        GM=sum(su[i].*(Z[:,i]*(Z[:,i]')-I) for i in 1:n)
        Gσ=tr(GM)/d
        GB=GM-Gσ*I
        μ=μ.+ημ*σ*B*Gδ
        σ=σ*exp(ησ/2*Gσ)
        B=B*exp(ηB/2 .* GB)
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
