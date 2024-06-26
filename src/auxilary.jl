######################
# Auxilary functions #
######################
"""
    dist(i::Int, j::Int, embed::Array{Float64,2})

Calculates Euclidian distance between two vectors from embedding array.

**Arguments**
* `v1::Int` index of first vector
* `v2::Int` index of second vector
* `embed::Array{Float64,2}` graph embedding array
"""
function dist(v1::Int, v2::Int, embed::Array{Float64,2})
    if v1 == v2
        return 0.0
    else
        return @inbounds sqrt(sum(i -> (embed[v1, i] - embed[v2, i])^2, axes(embed, 2)))
    end
end

"""
    JS(vC::Vector{Float64}, vB::Vector{Float64},
        vI::Vector{Int}, internal::Int, vLen::Int)

Jensen-Shannon divergence with Dirichlet-like prior.

**Arguments**
* `vC::Vector{Float64}` first distribution of edges within and between communities
* `vB::Vector{Float64}` second distribution of edges within and between communities
* `vI::Vector{Int}` indicator of internal (1) and external (0) edges w.r.t. communities, if empty compute overall JS distance
* `internal::Int` internal JS distance switch, if 1 return internal, else return external
"""
function JS(vC::Vector{Float64}, vB::Vector{Float64},
    vI::AbstractVector{Bool}, internal::Bool)

    if !isempty(vI)
        inter = vI .== internal
        vect_p = vC[inter]
        vect_q = vB[inter]
    else
        vect_p = copy(vC)
        vect_q = copy(vB)
    end
    sp_1 = sum(vect_p) + length(vect_p)
    sp_2 = sum(vect_q) + length(vect_q)
    @. vect_p = (vect_p + 1) / sp_1
    @. vect_q = (vect_q + 1) / sp_2
    vect_m = (vect_p .+ vect_q) ./ 2
    f = sum(@. vect_p * log(vect_p / vect_m) + vect_q * log(vect_q / vect_m))
    return f / 2
end

"""
Change 2-dimensional into 1-dim index
"""
function idx(n::Int, i::Int, j::Int)
    return n * (i - 1) - (i - 1) * (i - 2) ÷ 2 + j - i + 1
end

parse_flag(flag::String) = !isnothing(findfirst(==(flag), ARGS)) ? true : false

function parseargs()
    methods = Dict("rss" => split_cluster_rss,
        "rss2" => split_cluster_rss2,
        "size" => split_cluster_size,
        "diameter" => split_cluster_diameter)
    try
        ###Flags###

        # Check if calculations should be verbose
        verbose = parse_flag("-v")
        # Check if provided graph is directed
        directed = parse_flag("-d")
        # Check if global score should be done with split
        split = parse_flag("--split-global")

        ###Edgelist###

        idx = findfirst(==("-g"), ARGS)
        @assert !isnothing(idx) "Edgelist file is required"
        fn_edges = ARGS[idx+1]
        @assert isfile(fn_edges) "$fn_edges is not a file"

        # Read edges
        edges = readdlm(fn_edges, Float64)
        rows, no_cols = size(edges)
        verbose && println("$no_cols columns and $rows rows in edgelist file.")

        # Validate file structure
        @assert no_cols == 2 || no_cols == 3 "Expected 2 or 3 columns in edgelist file"
        v_min = minimum(edges[:, 1:2])
        @assert v_min == 0 || v_min == 1 "Vertices should be either 0-based or 1-based"

        # Make vertices 1-based
        if v_min == 0.0
            edges[:, 1:2] .+= 1.0
        end
        no_vertices = Int(maximum(edges[:, 1:2]))
        verbose && println("Graph contains $no_vertices vertices")

        # If graph is unweighted, add unit weights
        # Compute vertices weights
        vweight = zeros(no_vertices)
        eweights = no_cols == 2 ? ones(rows) : edges[:, 3]
        edges = convert.(Int, edges[:, 1:2])
        for i in 1:rows
            vweight[edges[i, 1]] += eweights[i]
            vweight[edges[i, 2]] += eweights[i]
        end
        verbose && println("Done preparing edgelist and vertices weights")

        ###Communities###

        idx = findfirst(==("-c"), ARGS)
        if !isnothing(idx)
            fn_comm = ARGS[idx+1]
        else
            no_cols == 2 ? louvain_clust(v_min, fn_edges) : louvain_clust(fn_edges, edges, eweights)
            fn_comm = fn_edges * ".ecg"
        end
        comm = readdlm(fn_comm, Int)
        comm_rows, no_cols = size(comm)

        # Validate file structure
        @assert comm_rows == no_vertices "No. communities ($comm_rows) differ from no. nodes ($no_vertices)"
        @assert no_cols == 1 || no_cols == 2 "Expected 1 or 2 columns in communities file, but encountered $no_cols."
        # 2 columns file - sort by first column and extract only second column
        if no_cols == 2
            comm = comm[sortperm(comm[:, 1]), 2]
            comm = reshape(comm, size(comm)[1], 1)
        end
        c_min = minimum(comm)
        @assert c_min == 0 || c_min == 1 "Communities should be either 0-based or 1-based, but are $c_min based."

        # Make communities 1-based
        if c_min == 0
            comm .+= 1
        end
        verbose && println("Done preparing communities.")

        ###Embedding###

        idx = findfirst(==("-e"), ARGS)
        @assert !isnothing(idx) "Embedding file is required"
        fn_embed = ARGS[idx+1]
        @assert isfile(fn_embed) "$fn_embed is not a file"

        # Read embedding
        embedding = []
        try
            embedding = readdlm(fn_embed, Float64)
        catch
            verbose && println("Embedding in node2vec format. Loading without first line.")
            embedding = readdlm(fn_embed, Float64, skipstart=1)
        end
        # Validate file
        @assert no_vertices == size(embedding, 1) "No. rows in embedding and no. vertices in a graph differ."

        # If embedding contains indices in first column, sort by it and remove the column
        try
            order = convert.(Int, embedding[:, 1])
            verbose && println("Sorting embedding by first column")
            embedding = embedding[sortperm(order), 2:end]
        catch

        end
        verbose && println("Done preparing embedding.")

        ###Landmarks###

        # Transform communities
        clusters = Dict{Any,Vector{Int}}()

        landmarks = -1
        idx = findfirst(==("-l"), ARGS)
        if !isnothing(idx)
            try
                landmarks = parse(Int, ARGS[idx+1])
            catch
                landmarks = round(Int, 4 * sqrt(no_vertices))
                @info "Using $landmarks landmarks"
            end
        end

        idx = findfirst(==("-f"), ARGS)
        if !isnothing(idx)
            forced = parse(Int, ARGS[idx+1])
            landmarks = landmarks == -1 ? 1 : landmarks
        else
            forced = 4
        end

        if no_vertices >= 10000 && isnothing(findfirst(==("--force-exact"), ARGS)) && landmarks == -1
            landmarks = max(round(Int, 4 * sqrt(no_vertices)), 4 * maximum(comm))
            @info "Number of vertices is equal or higher than 10 000. Automatically switching to approximate algortihm with $landmarks landmarks. If you want to force exact algorithm use --force-exact flag."
        end
        # Provide clusters membership only if generating landmarks
        if landmarks != -1
            for (i, c) in enumerate(comm[:, 1])
                if haskey(clusters, c)
                    push!(clusters[c], i)
                else
                    clusters[c] = [i]
                end
            end
            clusters = collect(values(clusters))
        end

        idx = findfirst(==("--seed"), ARGS)
        seed = !isnothing(idx) ? parse(Int, ARGS[idx+1]) : -1

        idx = findfirst(==("--samples-local"), ARGS)
        samples = !isnothing(idx) ? parse(Int, ARGS[idx+1]) : 10000

        idx = findfirst(==("-m"), ARGS)
        method_str = !isnothing(idx) ? lowercase(strip(ARGS[idx+1])) : "rss"

        method = methods[method_str]
        return edges, eweights, vweight, comm, clusters, embedding, verbose, landmarks, forced, method, directed, split, seed, samples
    catch e
        showerror(stderr, e)
        println("\n\nUsage:")
        println("\tjulia CGE_CLI.jl -g edgelist -e embedding [-c communities] [--seed seed] [--samples-local samples] [-v] [-d] [--split-global] [-l [landmarks]] [-f [forced]] [--force-exact] [-m method]")
        println("\nParameters:")
        println("edgelist: rows should contain two whitespace separated vertices ids (edge) and optional weights in third column")
        println("embedding: rows should contain whitespace separated embeddings of vertices")
        println("communities: rows should contain cluster identifiers of vertices with optional vertices ids in the first column")
        println("if no file is given communities are calculated with Louvain algorithm")
        println("seed: RNG seed for local measure sampling")
        println("samples: no. samples to draw for local score calculation")
        println("-v: flag for debugging messages")
        println("-d: flag for usage of directed framework")
        println("--split-global: flag for using splitted global score; kept for backward compatibility")
        println("landmarks: required number of landmarks; 4*sqrt(no.vertices) by default")
        println("forced: required number of forced splits of a cluster; 4 by default")
        println("if both 'landmarks' and 'multiplier' are provided the higher value is taken")
        println("method: one of:")
        println("\t* rss:      minimize maximum residual sum of squares when doing a cluster split")
        println("\t* rss2:     minimize maximum residual sum of squares when doing a cluster split (slower)")
        println("\t* size:     make clusters have approximately the same size after a cluster split")
        println("\t* diameter: make clusters have approximately the same diameter along first " *
                "principal component after a cluster split")
        println("(note that always a cluster with largest residual sum of squares is selected for splitting)")
        exit(1)
    end
end
