function randchi(rn::AbstractRNG,::Type{T},ν::Integer) where T
    ν < 0 && throw(ArgumentError("ν (DOF) must be atleast 0."))
    s=zero(T)
    for i in 1:ν
        s+=abs2(randn(rn,T))
    end
    s
end

function separable_nes(f,x0::AbstractVector{T},σ::AbstractVector{T},lr::ParamTuple{T}) where T
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
            else
                rescaling=sqrt(randchi(rn,T,N))/norm(ϵ[j])
                ϵ[j].*=rescaling
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
        x .+= lr.ημ .* σ .* ∇f
        σ .*= exp.(lr.ησ ./ 2 .*∇fσ)
    end
    (sol=x, cost=f(x))
end

function separable_nes(f,μ::AbstractVector{T},σ::T,lr::ParamTuple{T}) where T
    separable_nes(f,μ,fill(σ,length(μ)),lr)
end

function separable_nes(f,μ::AbstractVector{T},σ) where T
    separable_nes(f,μ,σ,LearningRates(separable_nes,T(length(μ))))
end

function LearningRates(::typeof(separable_nes),d)
    ημ=d/d
    ησ=3/5*(3+log(d))/(d*√d)
    (ημ=ημ,ησ=ησ)
end
