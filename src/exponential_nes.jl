struct xNES{T}
    ημ::T
    ησ::T
    ηB::T
    σtol::T
    samples::Int
    function xNES{T}(d::Integer;P...) where T
        ημ=T(1)
        ηB=ησ=T( (9+3*log(d))/(5*d*√d) )
        samples=4 + ceil(Int, log(3*d))
        σtol=T(1e-8)
        haskey(P,:ημ) && (ημ=P[:ημ])
        haskey(P,:ησ) && (ησ=P[:ησ])
        haskey(P,:ηB) && (ησ=P[:ηB])
        haskey(P,:samples) && (samples=P[:samples])
        haskey(P,:σtol) && (σtol=P[:σtol])
        new{T}(ημ,ησ,ηB,σtol,samples)
    end
end

function δ(T,i,j)
    ifelse(i==j,one(T),zero(T))
end
function exponential_nes(f,μ0::AbstractVector{T},A::AbstractMatrix{T},params::xNES{T}) where T
    n=params.samples
    ημ=params.ημ
    ησ=params.ησ
    ηB=params.ηB
    σtol=params.σtol
    μ=copy(μ0)
    d=length(μ)
    Z=Array{T}(undef,d,n)
    σ=abs(det(A))^(1/d)
    F=fill(f(μ),n)
    B=A./σ
    idx=collect(1:n)
    u=utility_function{T}(n)
    Gδ=zeros(T,d)
    GM=zeros(T,d,d)
    while σ>σtol
        randn!(Z)
        for i in 1:n
            F[i]=f(μ .+ σ .* B * Z[:,i])
        end
        sort!(idx,by=i->F[i])
        #Gδ=sum(u[i] .* Z[:,idx[i]] for i in 1:n)
        #GM=sum(u[i] .* (Z[:,idx[i]] * (Z[:,idx[i]]') - I) for i in 1:n)
        Gδ.=0
        GM.=0
        for i in 1:n
            j=idx[i]
            for k in 1:d
                Gδ[k]+=u[i] * Z[k,j]
            end
            for k1 in 1:d, k2 in 1:d
                GM[k1,k2]+=u[i] * (Z[k1,j] * Z[k2,j] - δ(T,k1,k2))
            end
        end
        Gσ=tr(GM) / d
        #=
        GB=GM - Gσ * I
        μ=μ .+ ημ * σ *B *Gδ
        σ=σ * exp(ησ/2 * Gσ)
        B=B * exp(ηB/2 .* GB)=#
        Gσ=tr(GM) / d
        for i in 1:d
            GM[i,i]-=Gσ
        end
        μ=μ .+ ημ * σ *B *Gδ
        σ=σ * exp(ησ/2 * Gσ)
        B=B * exp(ηB/2 .* GM)
    end
    return (sol=μ, cost=f(μ))
end

function exponential_nes(f,μ::AbstractVector{T},σ::AbstractVector{T},params::xNES{T}) where T
    A=diagm(σ)
    exponential_nes(f,μ,A,params)
end

function exponential_nes(f,μ::AbstractVector{T},σ::T,params::xNES{T}) where T
    exponential_nes(f,μ,fill(σ,length(μ)),params)
end

function exponential_nes(f,μ::AbstractVector{T},σ; P...) where T
    exponential_nes(f,μ,σ,xNES{T}(length(μ);P...))
end
