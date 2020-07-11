using CGE

ARGS=["-g","example/100k.edgelist","-c","example/100k.ecg","-e","example/100k.embedding","-o","testfile","-l","200","-f","1","-m","rss"];
edges, weights, vweights, comm, clusters, embed, asis, verbose, land, forced, method = parseargs(ARGS)
distances = zeros(length(vweights))
if land != -1
    distances, embed, comm, edges, weights  = landmarks(edges, weights, vweights, clusters, comm,
                                                        embed, verbose, land, forced, method)
end
results = wGCL(edges, weights, comm, embed, distances, verbose)
println(results)