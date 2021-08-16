using CGE

edges, weights, vweights, comm, clusters, embed, verbose, land, forced, method, directed = parseargs()
distances = zeros(length(vweights))
init_edges = Array{Int,2}(undef,0,0)
init_embed = Array{Float64,2}(undef,0,0)
lweight = Float64[]
v_to_l = Int[]
if land != -1
    init_edges = copy(edges)
    init_embed = copy(embed)
    distances, embed, comm, edges, weights, lweight, v_to_l = landmarks(edges, 
        weights, vweights, clusters, comm, embed, verbose, land, forced, method)
end
if directed
    results = wGCL_directed(edges, weights, comm, embed, distances, vweights, lweight, v_to_l, init_edges, init_embed, verbose)
else
    results = wGCL(edges, weights, comm, embed, distances, vweights, lweight, v_to_l, init_edges, init_embed, verbose)
end
println(results)