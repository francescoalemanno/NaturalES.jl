mutable struct sNES_state{T}
    ημ::T
    ησ::T
    samples::Int
    termination::Termination{T}
    function sNES_state{T}(d::Integer;P...) where T
        ημ=T(1)
        ησ=T( (3+log(d))/(5*√d) )
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

struct sNES <: OptMethod end

function separable_nes(f,x0::AbstractVector{T},σ::AbstractVector{T},params::sNES_state{T}) where T
    samples=params.samples
    samples>3 || error("atleast 3 samples are required for sNES to work.")
    ημ = params.ημ
    ησ = params.ησ
    N=length(x0)
    length(σ) == N || error("the length of 'σ' must be equal to the length of 'x'.")
    x = copy(x0)
    Z = zero(T)
    O = one(T)
    F = fill(f(x0),samples)
    ∇f = fill(Z,N)
    ∇fσ = fill(Z,N)
    ϵ = map(i->zeros(T,N), 1:samples)
    u = utility_function{T}(samples)
    idx = collect(1:samples)
    tmp_x = copy(x0)
    for i in 1:samples
        randn!(ϵ[i])
        tmp_x .= x .+ σ .* ϵ[i]
        F[i] = f(tmp_x)
    end
    sortperm!(idx,F)
    iterations=0
    while !termination_criteria(params.termination,geo_mean(σ),iterations)
        for i in 1:samples
            j = idx[i]
            if i >= (samples ÷ 2)
                randn!(ϵ[j])
            else
                rescaling = sqrt(randchisq(T, N)) / norm(ϵ[j])
                ϵ[j] .*= rescaling
            end
            tmp_x .= x .+ σ .* ϵ[j]
            F[j] = f(tmp_x)
        end
        sortperm!(idx,F)
        ∇f .= Z
        ∇fσ .= Z
        for i in 1:samples
            j = idx[i]
            ∇f .+= u[i] .* ϵ[j]
            ∇fσ .+= u[i] .* (ϵ[j] .* ϵ[j] .- O)
        end
        x .+= ημ .* σ .* ∇f
        σ .*= exp.(ησ ./ 2 .* ∇fσ)
        iterations+=1
    end
    (sol = x, cost = f(x))
end

function separable_nes(f,μ::AbstractVector{T},σ::T,params::sNES_state{T}) where T
    separable_nes(f,μ,fill(σ,length(μ)),params)
end

function separable_nes(f,μ::AbstractVector{T},σ; P...) where T
    separable_nes(f,μ,σ,sNES_state{T}(length(μ);P...))
end

pickmethod(::Type{sNES}) = separable_nes
