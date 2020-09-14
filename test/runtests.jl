using Test
using CGE

for e in ["-g","test.edgelist","-c","test.ecg","-e","test.embedding","-l","20","-f","1","-m","rss"]
    push!(Base.ARGS,e)
end
edges, weights, vweights, comm, clusters, embed, asis, verbose, land, forced, method = parseargs()

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

distances, lembed, lcomm, ledges, lweights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, method);

@testset "landmarks-rss" begin

    @test typeof(ledges) == Array{Int,2}
    @test minimum(ledges) == 1
    @test typeof(lweights) == Array{Float64,1}
    @test typeof(lcomm) == Array{Int,2}
    @test size(lcomm,2) == 1
    @test typeof(distances) == Array{Float64,1}
end

distances, lembed, lcomm, ledges, lweights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, CGE.split_cluster_rss2);

@testset "landmarks-rss2" begin

@test typeof(ledges) == Array{Int,2}
@test minimum(ledges) == 1
@test typeof(lweights) == Array{Float64,1}
@test typeof(lcomm) == Array{Int,2}
@test size(lcomm,2) == 1
@test typeof(distances) == Array{Float64,1}
end

distances, lembed, lcomm, ledges, lweights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, CGE.split_cluster_size);

@testset "landmarks-size" begin

@test typeof(ledges) == Array{Int,2}
@test minimum(ledges) == 1
@test typeof(lweights) == Array{Float64,1}
@test typeof(lcomm) == Array{Int,2}
@test size(lcomm,2) == 1
@test typeof(distances) == Array{Float64,1}
end

distances, lembed, lcomm, ledges, lweights  = landmarks(edges, weights, vweights, clusters, comm,
embed, verbose, land, forced, CGE.split_cluster_diameter);

@testset "landmarks-diameter" begin

@test typeof(ledges) == Array{Int,2}
@test minimum(ledges) == 1
@test typeof(lweights) == Array{Float64,1}
@test typeof(lcomm) == Array{Int,2}
@test size(lcomm,2) == 1
@test typeof(distances) == Array{Float64,1}
end

results = wGCL(ledges, lweights, lcomm, lembed, distances, verbose);

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