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
    if v1==v2
        return 0.0
    else
        return @inbounds sqrt(sum(i -> (embed[v1, i] - embed[v2, i]) ^ 2, axes(embed, 2)))
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
    f = sum(@. vect_p*log(vect_p/vect_m)+vect_q*log(vect_q/vect_m))
    return f / 2
end

"""
Change 2-dimensional into 1-dim index
"""
function idx(n::Int, i::Int, j::Int)
  return n*(i-1) - (i-1)*(i-2) ÷ 2 +j-i+1
end

function parseargs()
    methods = Dict("rss" => split_cluster_rss,
                   "rss2" => split_cluster_rss2,
                   "size" => split_cluster_size,
                   "diameter" => split_cluster_diameter)
    try
        # Optional arguments
        ##Flags
        asis = !isnothing(findfirst(==("-a"),ARGS)) ? true : false
        # Check if calculations should be verbose
        verbose = !isnothing(findfirst(==("-v"),ARGS)) ? true : false

        # Check for required arguments: -g graph_edgelist -c communities -e embedding -o outfile
        @assert length(ARGS) >= 6

        # Load obligatory files
        ################
        ## Graph edges #
        ################

        idx = findfirst(==("-g"),ARGS)
        @assert !isnothing(idx) "Edges list file is required"
        fn_edges = ARGS[idx+1]

        # read edges
        edges = readdlm(fn_edges, Float64)
        rows, no_cols = size(edges)
        verbose && println("$no_cols columns in graph edgelist file")

        # Validate file structure
        @assert no_cols==2 || no_cols==3 "Expected 2 or 3 columns"
        v_min = minimum(edges[:,1:2])
        @assert v_min==0 || v_min==1 "Vertices should be either 0-based or 1-based"

        # make vertices 1-based
        if v_min == 0.0
            edges[:,1:2] .+= 1.0
        end
        no_vertices = Int(maximum(edges[:,1:2]))
        verbose && println("Vertices from $v_min to $no_vertices")

        # if unweighted, add unit weights
        # compute vertex weights
        if no_cols == 2
            edges = convert.(Int,edges)
            eweights = ones(rows)
            vweight = zeros(no_vertices)
            for i in 1:rows
                vweight[edges[i,1]]+=1.0
                vweight[edges[i,2]]+=1.0
            end
        else
            edges = convert.(Int,edges[:,1:2])
            eweights = edges[:,3]
            vweight = zeros(no_vertices)
            for i in 1:rows
                vweight[edges[i,1]]+=eweights[i]
                vweight[edges[i,2]]+=eweights[i]
            end
        end
        verbose && println("Done preparing edgelist and vertex weight")

        ################
        ## Communities #
        ################

        idx = findfirst(==("-c"),ARGS)
        @assert !isnothing(idx) "Communities file is required"
        fn_comm = ARGS[idx+1]

        # Read communities
        comm = readdlm(fn_comm,Int)
        comm_rows, no_cols = size(comm)

        # Validate file structure
        @assert no_cols==1 "Expected 1 column"
        v_min = minimum(comm[:,1])
        @assert v_min==0 || v_min==1 "Communities should be either 0-based or 1-based"
        c_min=minimum(comm)

        # make communities 1-based
        if c_min == 0
            comm[:,1] .+=1
        end
        verbose && println("Done preparing communities")

        ##############
        ## Embedding #
        ##############

        idx = findfirst(==("-e"),ARGS)
        @assert !isnothing(idx) "Embedding file is required"
        fn_embed = ARGS[idx+1]

        # Read embedding
        embedding = readdlm(fn_embed,Float64)

        # Validate file
        @assert comm_rows == size(embedding, 1) "No. rows in embedding and communities files differ"

        # if embedding contains index in first column, sort by it and remove column
        if asis == 0
            embedding = embedding[sortperm(embedding[:,1]),2:end]
        end
        verbose && println("done preparing embedding")

        #####################
        ## Output file name #
        #####################
        # idx = findfirst(==("-o"),ARGS)
        # @assert !isnothing(idx)
        # outfile = ARGS[idx+1]

        #############
        ##Landmarks #
        #############

        # Transform communities
        clusters = Dict{Any, Vector{Int}}()
        for (i, c) in enumerate(comm[:,1])
            if haskey(clusters, c)
                push!(clusters[c], i)
            else
                clusters[c] = [i]
            end
        end
        clusters = collect(values(clusters))

        idx = findfirst(==("-l"),ARGS)
        landmarks = !isnothing(idx) ? parse(Int, ARGS[idx+1]) : -1
        idx = findfirst(==("-f"),ARGS)
        forced = !isnothing(idx) ? parse(Int, ARGS[idx+1]) : -1
        idx = findfirst(==("-m"),ARGS)
        method_str = !isnothing(idx) ? lowercase(strip(ARGS[idx+1])) : "rss"
        method = methods[method_str]
        return edges, eweights, vweight, comm, clusters, embedding, asis, verbose, landmarks, forced, method
    catch e
        showerror(stderr, e)
        println("\n\nUsage:")
        println("\tjulia CGE.jl -g graph_edgelist -c communities -e embedding [-a -v] [-l landmarks -f forced -m method]")
        println("\nParameters:")
        println("graph_edgelist: rows should contain two vertices ids (edge) and optional weights")
        println("communities: rows should contain cluster identifiers of consecutive vertices")
        println("embedding: rows should contain whitespace separated locations of vertices in embedding")
        println("-a: flag for sorting embedding")
        println("-v: flag for debugging messages")
        println("landmarks: required number of landmarks")
        println("forced: required maximum number of forced splits of a cluster")
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

function parseargs(ARGS::Vector{String})
    methods = Dict("rss" => split_cluster_rss,
                   "rss2" => split_cluster_rss2,
                   "size" => split_cluster_size,
                   "diameter" => split_cluster_diameter)
    try
        # Optional arguments
        ##Flags
        asis = !isnothing(findfirst(==("-a"),ARGS)) ? true : false
        # Check if calculations should be verbose
        verbose = !isnothing(findfirst(==("-v"),ARGS)) ? true : false

        # Check for required arguments: -g graph_edgelist -c communities -e embedding
        @assert length(ARGS) >= 6

        # Load obligatory files
        ################
        ## Graph edges #
        ################

        idx = findfirst(==("-g"),ARGS)
        @assert !isnothing(idx)
        fn_edges = ARGS[idx+1]

        # read edges
        edges = readdlm(fn_edges, Float64)
        rows, no_cols = size(edges)
        verbose && println("$no_cols columns in graph edgelist file")

        # Validate file structure
        @assert no_cols==2 || no_cols==3 "Expected 2 or 3 columns"
        v_min = minimum(edges[:,1:2])
        @assert v_min==0 || v_min==1 "Vertices should be either 0-based or 1-based"

        # make vertices 1-based
        if v_min == 0.0
            edges[:,1:2] .+= 1.0
        end
        no_vertices = Int(maximum(edges[:,1:2]))
        verbose && println("Vertices from $v_min to $no_vertices")

        # if unweighted, add unit weights
        # compute vertex weights
        if no_cols == 2
            edges = convert.(Int,edges)
            eweights = ones(rows)
            vweight = zeros(no_vertices)
            for i in 1:rows
                vweight[edges[i,1]]+=1.0
                vweight[edges[i,2]]+=1.0
            end
        else
            edges = convert.(Int,edges[:,1:2])
            eweights = edges[:,3]
            vweight = zeros(no_vertices)
            for i in 1:rows
                vweight[edges[i,1]]+=eweights[i]
                vweight[edges[i,2]]+=eweights[i]
            end
        end
        verbose && println("Done preparing edgelist and vertex weight")

        ################
        ## Communities #
        ################

        idx = findfirst(==("-c"),ARGS)
        @assert !isnothing(idx)
        fn_comm = ARGS[idx+1]

        # Read communities
        comm = readdlm(fn_comm,Int)
        comm_rows, no_cols = size(comm)

        # Validate file structure
        @assert no_cols==1 "Expected 1 column"
        v_min = minimum(comm[:,1])
        @assert v_min==0 || v_min==1 "Communities should be either 0-based or 1-based"
        c_min=minimum(comm)

        # make communities 1-based
        if c_min == 0
            comm[:,1] .+=1
        end
        verbose && println("Done preparing communities")

        ##############
        ## Embedding #
        ##############

        idx = findfirst(==("-e"),ARGS)
        @assert !isnothing(idx)
        fn_embed = ARGS[idx+1]

        # Read embedding
        embedding = readdlm(fn_embed,Float64)

        # Validate file
        @assert comm_rows == size(embedding, 1) "No. rows in embedding and communities files differ"

        # if embedding contains index in first column, sort by it and remove column
        if asis == 0
            embedding = embedding[sortperm(embedding[:,1]),2:end]
        end
        verbose && println("done preparing embedding")

        #####################
        ## Output file name #
        #####################
        # idx = findfirst(==("-o"),ARGS)
        # @assert !isnothing(idx)
        # outfile = ARGS[idx+1]

        #############
        ##Landmarks #
        #############

        # Transform communities
        clusters = Dict{Any, Vector{Int}}()
        for (i, c) in enumerate(comm[:,1])
            if haskey(clusters, c)
                push!(clusters[c], i)
            else
                clusters[c] = [i]
            end
        end
        clusters = collect(values(clusters))

        idx = findfirst(==("-l"),ARGS)
        landmarks = !isnothing(idx) ? parse(Int, ARGS[idx+1]) : -1
        idx = findfirst(==("-f"),ARGS)
        forced = !isnothing(idx) ? parse(Int, ARGS[idx+1]) : -1
        idx = findfirst(==("-m"),ARGS)
        method_str = !isnothing(idx) ? lowercase(strip(ARGS[idx+1])) : "rss"
        method = methods[method_str]
        return edges, eweights, vweight, comm, clusters, embedding, asis, verbose, landmarks, forced, method
    catch e
        showerror(stderr, e)
        println("\n\nUsage:")
        println("\tjulia CGE.jl -g graph_edgelist -c communities -e embedding [-a -v] [-l landmarks -f forced -m method]")
        println("\nParameters:")
        println("graph_edgelist: rows should contain two vertices ids (edge) and optional weights")
        println("communities: rows should contain cluster identifiers of consecutive vertices")
        println("embedding: rows should contain whitespace separated locations of vertices in embedding")
        println("-a: flag for sorting embedding")
        println("-v: flag for debugging messages")
        println("landmarks: required number of landmarks")
        println("forced: required maximum number of forced splits of a cluster")
        println("method: one of:")
        println("\t* rss:      minimize maximum residual sum of squares when doing a cluster split")
        println("\t* rss2:     minimize maximum residual sum of squares when doing a cluster split (slower)")
        println("\t* size:     make clusters have approximately the same size after a cluster split")
        println("\t* diameter: make clusters have approximately the same diameter along first " *
                "principal component after a cluster split")
        println("(note that always a cluster with largest residual sum of squares is selected for splitting)")
    end
end