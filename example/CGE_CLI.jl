using CGE

edges, weights, vweights, comm, clusters, embed, verbose, land, forced, method, directed = parseargs()
distances = zeros(length(vweights))
if land != -1
    distances, embed, comm, edges, weights  = landmarks(edges, weights, vweights, clusters, comm,
                                                        embed, verbose, land, forced, method)
end
println("distances len:",length(distances))
println("embed len:",length(embed))
println("comm len:",length(comm))
println("edges len:",length(edges))
println("weights len:",length(weights))
if directed
    results = wGCL_directed(edges, weights, comm, embed, distances, verbose)
else
    results = wGCL(edges, weights, comm, embed, distances, verbose)
end
println(results)