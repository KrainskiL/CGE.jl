module CGE

using DelimitedFiles
using Statistics
using LinearAlgebra
using louvain_jll

#auxilary
export parseargs

#landmarks
export landmarks

#divergence
export wGCL
export wGCL_directed

#clustering
export louvain_clust

# Include package code
include("auxilary.jl")
include("landmarks.jl")
include("divergence.jl")
include("clustering.jl")
end
