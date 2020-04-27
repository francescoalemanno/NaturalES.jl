mutable struct xNES_state{T}
    ημ::T
    ησ::T
    samples::Int
    termination::Termination{T}
    function xNES_state{T}(d::Integer;P...) where T
        ημ=T(1)
        ησ=T( (9+3*log(d))/(5*d*√d) )
        samples=4 + ceil(Int, log(3*d))
        Term=Termination{T}()
        S=new{T}(ημ,ησ,samples,Term)
        for k in keys(P)
            hasfield(sNES_state{T},k) && setfield!(S,k,P[k])
            hasfield(Termination{T},k) && setfield!(Term,k,P[k])
        end
        S
    end
end

struct xNES <: OptMethod end

function δ(T,i,j)
    ifelse(i==j,one(T),zero(T))
end

function exponential_nes(f,μ0::AbstractVector{T},A::AbstractMatrix{T},params::xNES_state{T}) where T
    n=params.samples
    ημ=params.ημ
    ηB=ησ=params.ησ
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
    tmp_μ=copy(μ)
    iterations=0
    @inbounds while !termination_criteria(params.termination,σ,iterations)
        randn!(Z)
        for i in 1:n
            #F[i] = f(μ .+ σ .* B * Z[:,i])
            mul!(tmp_μ,B,view(Z,:,i))
            tmp_μ .= μ .+ tmp_μ .* σ
            F[i] = f(tmp_μ)
        end
        sortperm!(idx,F)
        #Gδ=sum(u[i] .* Z[:,idx[i]] for i in 1:n)
        #GM=sum(u[i] .* (Z[:,idx[i]] * (Z[:,idx[i]]') - I) for i in 1:n)
        let i=1
            j = idx[i]
            for k2 in 1:d
                Gδ[k2] = u[i] * Z[k2,j]
                for k1 in 1:d
                    GM[k1,k2] = u[i] * (Z[k1,j] * Z[k2,j] - δ(T,k1,k2))
                end
            end
        end
        for i in 2:n
            j = idx[i]
            for k2 in 1:d
                Gδ[k2] += u[i] * Z[k2,j]
                for k1 in 1:d
                    GM[k1,k2] += u[i] * (Z[k1,j] * Z[k2,j] - δ(T,k1,k2))
                end
            end
        end
        #=
        Gσ = tr(GM) / d
        GB=GM - Gσ * I
        μ=μ .+ ημ * σ *B *Gδ
        σ=σ * exp(ησ/2 * Gσ)
        B=B * exp(ηB/2 .* GB)
        =#
        Gσ = tr(GM) / d
        for i in 1:d
            GM[i,i]-=Gσ
        end
        GM .*= ηB/2
        Gσ *= ησ/2
        mul!(tmp_μ,B,Gδ)
        μ .+= ημ .* σ .* tmp_μ
        σ *= exp(Gσ)
        mul!(GM,B,exp(GM))
        B.=GM
        iterations+=1
    end
    return (sol=μ, cost=f(μ))
end

function exponential_nes(f,μ::AbstractVector{T},σ::AbstractVector{T},params::xNES_state{T}) where T
    A=diagm(σ)
    exponential_nes(f,μ,A,params)
end

function exponential_nes(f,μ::AbstractVector{T},σ::T,params::xNES_state{T}) where T
    exponential_nes(f,μ,fill(σ,length(μ)),params)
end

function exponential_nes(f,μ::AbstractVector{T},σ; P...) where T
    exponential_nes(f,μ,σ,xNES_state{T}(length(μ);P...))
end

pickmethod(::Type{xNES}) = exponential_nes
