using Test
using CGE

for e in ["-g","test.edgelist","-c","test1col.ecg","-e","test_n2v.embedding","-l","20","-f","1","-m","rss"]
    push!(ARGS,e)
end
edges, weights, vweights, comm, clusters, embed, verbose, land, forced, method = parseargs()

empty!(ARGS)
for e in ["-g","test.edgelist","-c","test2col.ecg","-e","test_ordered.embedding","-l","20","-f","1","-m","rss"]
    push!(ARGS,e)
end
edges, weights, vweights, comm, clusters, embed, verbose, land, forced, method = parseargs()

empty!(ARGS)
for e in ["-g","test_weights.edgelist","-c","test2col.ecg","-e","test_unordered.embedding","-l","20","-f","1","-m","rss"]
    push!(ARGS,e)
end
edges, weights, vweights, comm, clusters, embed, verbose, land, forced, method = parseargs()

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

@test typeof(verbose) == Bool

@test typeof(land) == Int
@test typeof(forced) == Int

end

distances, lembed, lcomm, ledges, lweights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, method, false);

@testset "landmarks-rss" begin

    @test typeof(ledges) == Array{Int,2}
    @test minimum(ledges) == 1
    @test typeof(lweights) == Array{Float64,1}
    @test typeof(lcomm) == Array{Int,2}
    @test size(lcomm,2) == 1
    @test typeof(distances) == Array{Float64,1}
end

distances, lembed, lcomm, ledges, lweights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, CGE.split_cluster_rss2,false);

@testset "landmarks-rss2" begin

@test typeof(ledges) == Array{Int,2}
@test minimum(ledges) == 1
@test typeof(lweights) == Array{Float64,1}
@test typeof(lcomm) == Array{Int,2}
@test size(lcomm,2) == 1
@test typeof(distances) == Array{Float64,1}
end

distances, lembed, lcomm, ledges, lweights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, CGE.split_cluster_size,false);

@testset "landmarks-size" begin

@test typeof(ledges) == Array{Int,2}
@test minimum(ledges) == 1
@test typeof(lweights) == Array{Float64,1}
@test typeof(lcomm) == Array{Int,2}
@test size(lcomm,2) == 1
@test typeof(distances) == Array{Float64,1}
end

distances, lembed, lcomm, ledges, lweights, lvweights, v_to_l  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, CGE.split_cluster_diameter,false);

@testset "landmarks-diameter" begin

@test typeof(ledges) == Array{Int,2}
@test minimum(ledges) == 1
@test typeof(lweights) == Array{Float64,1}
@test typeof(lcomm) == Array{Int,2}
@test size(lcomm,2) == 1
@test typeof(distances) == Array{Float64,1}
end

results = wGCL(ledges, lweights, lcomm, lembed, distances, lvweights,
    Float64[], Int64[], Matrix{Int}(undef,0,0), Float64[],Matrix{Float64}(undef,0,0),false, 42, 10000,verbose);

@testset "wgcl" begin

@test typeof(results) == Array{Float64,1}
@test results[1] <=10.0

end

@testset "clustering" begin

louvain_clust("test.edgelist")
@test isfile("test.edgelist.ecg")
isfile("test.edgelist.ecg") && rm("test.edgelist.ecg")

louvain_clust("testw.edgelist",edges, weights)
@test isfile("testw.edgelist.ecg")
isfile("testw.edgelist.ecg") && rm("testw.edgelist.ecg")
end
