module NaturalES
    using Random
    using LinearAlgebra

    const ParamTuple{T} = NamedTuple{S,NTuple{N,T}} where {S,N}

    include("separable_nes.jl")
    include("exponential_nes.jl")

    export separable_nes, exponential_nes
end # module
