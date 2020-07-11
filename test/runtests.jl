using Test
using CGE

argslist=["-g","test.edgelist","-c","test.ecg","-e","test.embedding","-l","20","-f","1","-m","rss"];
edges, weights, vweights, comm, clusters, embed, asis, verbose, land, forced, method = parseargs(argslist)

@testset "parsing" begin

@test typeof(edges) == Array{Int,2}
@test minimum(edges) == 1

@test typeof(weights) == Array{Float64,1}
@test typeof(vweights) == Array{Float64,1}

@test typeof(comm) == Array{Int,2}
@test size(comm,2) == 1
@test minimum(comm) == 1

@test typeof(clusters) == Vector{Vector{Int}}
@test typeof(embed) == Array{Float64,2}

@test typeof(asis) == Bool
@test typeof(verbose) == Bool

@test typeof(land) == Int
@test typeof(forced) == Int

end

distances, embed, comm, edges, weights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, method);

@testset "landmarks" begin

@test typeof(edges) == Array{Int,2}
@test minimum(edges) == 1

@test typeof(weights) == Array{Float64,1}

@test typeof(comm) == Array{Int,2}
@test size(comm,2) == 1

@test typeof(distances) == Array{Float64,1}
end

results = wGCL(edges, weights, comm, embed, distances, verbose);

@testset "wgcl" begin

@test typeof(results) == Array{Float64,1}
@test results[1] <=10.0

end
