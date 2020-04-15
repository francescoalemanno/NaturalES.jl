function exponential_nes(f,μ0::AbstractVector{T},A::AbstractMatrix{T},lr::ParamTuple{T}) where T
    μ=copy(μ0)
    d=length(μ)
    n=4+floor(Int,3*log(d))
    Z=Array{T}(undef,d,n)
    σ=abs(det(A))^(1/d)
    F=fill(f(μ),n)
    B=A./σ
    idx=collect(1:n)
    u=@. max(0,log(n/2+1)-log(idx))
    u=u./sum(u) .- 1/n

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
        μ=μ.+lr.ημ*σ*B*Gδ
        σ=σ*exp(lr.ησ/2*Gσ)
        B=B*exp(lr.ηB/2 .* GB)
    end
    return (sol=μ, cost=f(μ))
end

function exponential_nes(f,μ::AbstractVector{T},σ::AbstractVector{T},lr::ParamTuple{T}) where T
    A=diagm(σ)
    exponential_nes(f,μ,A,lr)
end

function exponential_nes(f,μ::AbstractVector{T},σ::T,lr::ParamTuple{T}) where T
    exponential_nes(f,μ,fill(σ,length(μ)),lr)
end

function exponential_nes(f,μ::AbstractVector{T},σ) where T
    exponential_nes(f,μ,σ,LearningRates(exponential_nes,T(length(μ))))
end

function LearningRates(::typeof(exponential_nes),d)
    ημ=d/d
    ησ=3/5*(3+log(d))/(d*√d)
    (ημ=ημ,ησ=ησ,ηB=ησ)
end
