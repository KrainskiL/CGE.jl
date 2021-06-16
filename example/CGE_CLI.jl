using CGE

edges, weights, vweights, comm, clusters, embed, verbose, land, forced, method = parseargs()
distances = zeros(length(vweights))
if land != -1
    distances, embed, comm, edges, weights  = landmarks(edges, weights, vweights, clusters, comm,
                                                        embed, verbose, land, forced, method)
end
results = wGCL(edges, weights, comm, embed, distances, verbose)
println(results)