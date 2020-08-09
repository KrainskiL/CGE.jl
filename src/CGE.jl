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

#clustering
export louvain

# Include package code
include("auxilary.jl")
include("landmarks.jl")
include("divergence.jl")
include("clustering.jl")
end
