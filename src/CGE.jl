module CGE

using DelimitedFiles
using Statistics
using LinearAlgebra

#auxilary
export parseargs

#landmarks
export landmarks

#divergence
export wGCL

# Include package code
include("auxilary.jl")
include("landmarks.jl")
include("divergence.jl")
end
