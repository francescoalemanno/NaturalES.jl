module NaturalES
using Random
using LinearAlgebra

include("utils.jl")
include("separable_nes.jl")
include("exponential_nes.jl")

export separable_nes, exponential_nes
end # module
