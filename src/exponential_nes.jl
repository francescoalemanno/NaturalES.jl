struct xNES{T}
    ημ::T
    ησ::T
    ηB::T
    σtol::T
    samples::Int
    function xNES{T}(d::Integer;P...) where T
        ημ=T(1)
        ηB=ησ=T(3/5*(3+log(d))/(d*√d))
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

    while σ>σtol
        randn!(Z)
        for i in 1:n
            F[i]=f(μ .+ σ .* B * Z[:,i])
        end
        su=u[invperm(sort(idx,by=i->F[i]))]
        Gδ=sum(su[i] .* Z[:,i] for i in 1:n)
        GM=sum(su[i] .* (Z[:,i] * (Z[:,i]') - I) for i in 1:n)
        Gσ=tr(GM) / d
        GB=GM - Gσ * I
        μ=μ .+ ημ * σ *B *Gδ
        σ=σ * exp(ησ/2 * Gσ)
        B=B * exp(ηB/2 .* GB)
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
