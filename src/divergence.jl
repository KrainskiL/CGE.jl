###############################
# Weighted Geometric Chung Lu #
###############################
"""
wGCL(edges::Array{Int,2}, weights::Vector{Float64}, comm::Array{Int,2},
            embed::Array{Float64,2}, distances::Vector{Float64}, verbose::Bool = false)

Calculates Weighted Geometric Chung-Lu model and divergence score for graph and embedding.

**Arguments**
* `edges::Array{Int,2}` array with edges definition (two whitespace separated vertices ids)
* `weights::Vector{Float64}` edges weights
* `comm::Array{Int,2}` assignment of vertices to communities
* `embed::Array{Float64,2}` array with vertices embeddings
* `distances::Vector{Float64}` distances between vertices
* `verbose::Bool` verbose switch, if true prints additional processing information
"""
function wGCL(edges::Array{Int,2}, weights::Vector{Float64}, comm::Matrix{Int},
            embed::Matrix{Float64}, distances::Vector{Float64}, verbose::Bool=false)
    # default values
    epsilon = 0.25
    delta = 0.001
    AlphaMax = 10.0
    AlphaStep = 0.25
    alpha_div_counter = 5

    no_vertices = maximum(edges)
    no_edges = size(edges,1)
    verbose && println("$no_vertices vertices and $no_edges edges")

    # Read communities
    @assert size(comm,1) == no_vertices
    n_parts = maximum(comm)

    # Compute degrees
    degree = zeros(no_vertices)
    for i in 1:no_edges
        a = weights[i]
        degree[edges[i,1]]+=a
        degree[edges[i,2]]+=a
    end

    # Compute C-vector
    vect_len = Int(n_parts*(n_parts+1)/2)
    vect_C = zeros(Float64, vect_len)
    vect_B = zeros(Float64, vect_len)

    for i in 1:no_edges
        a = weights[i]
        j,k = extrema([comm[edges[i,1]],comm[edges[i,2]]])
        vect_C[idx(n_parts, j, k)] += a
    end
    # Indicator - internal to a community
    vect_I = falses(vect_len)
    j=1
    for i in 1:n_parts
        vect_I[j] = true
        j+=(n_parts-i+1)
    end
    best_div = best_div_ext = best_div_int = typemax(Float64)
    best_alpha = -1.0
    dim = size(embed,2)
    verbose && println("Embedding has $dim dimensions")
    # Loop over alpha's

    # Compute Euclidean distance vector D[] given embed and alpha
    p_len = no_vertices*(no_vertices+1) ÷ 2
    D = zeros(Float64,p_len)
    @assert size(distances,1) == no_vertices
    for i in 1:no_vertices
        for j in i:no_vertices
            l = idx(no_vertices, i, j)
            if i == j
                D[l] = distances[i]
            else
                D[l] = dist(i,j,embed)
            end
        end
    end
    lo,hi = extrema(D)

    # Loop here - exclude Alpha = 0
    T = ones(no_vertices)
    for alpha in AlphaStep:AlphaStep:(AlphaMax+delta)
        # Apply kernel (g(dist))
        GD = zeros(Float64,p_len)
        for i in 1:no_vertices
            for j in i:no_vertices
                k = idx(no_vertices,i,j)
                GD[k] = (D[k]-lo)/(hi-lo) # normalize to [0,1]
                GD[k] = (1-GD[k])^alpha # transform w.r.t. alpha
            end
        end
        # Learn GCL model numerically
        diff = 1.0
        while diff > delta # stopping criterion
            S = zeros(no_vertices)
            for i in 1:no_vertices
                for j in i:no_vertices
                    tmp = T[i]*T[j]*GD[idx(no_vertices, i, j)]
                    S[i] += tmp
                    if i!=j S[j]+=tmp end
                end
            end
            f = 0.0
            for i in 1:no_vertices
                move = epsilon*T[i]*(degree[i]/S[i]-1.0)
                T[i]+=move
                f = max(f,abs(degree[i]-S[i])) # convergence w.r.t. degrees
            end
            diff = f
            verbose && println("diff = $diff")
        end
        # Compute probas P[]
        P = zeros(Float64,p_len)
        for i in 1:no_vertices
            for j in i:no_vertices
                P[k] = T[i]*T[j]*GD[idx(no_vertices,i,j)]
            end
        end

        # Compute B-vector given P[] and comm[]
        vect_B = zeros(vect_len)
        for i in 1:no_vertices
            for j in i:no_vertices
                k,l = extrema([comm[i],comm[j]])
                vect_B[Int(n_parts*(k-1)-(k-1)*(k-2)/2+l-k+1)] += P[idx(no_vertices,i,j)]
            end
        end
        x = JS(vect_C, vect_B, vect_I, true)
        y = JS(vect_C, vect_B, vect_I, false)
        f = (x+y)/2.0
        if f < best_div
            best_div = f
            best_alpha = alpha
            best_div_ext = x
            best_div_int = y
        else
            alpha_div_counter -= 1
            alpha_div_counter == 0 && break
        end
    end
    return [best_alpha, best_div, best_div_ext, best_div_int]
end

"""
wGCL_directed(edges::Array{Int,2}, weights::Vector{Float64}, comm::Array{Int,2},
            embed::Array{Float64,2}, distances::Vector{Float64}, verbose::Bool = false)

Calculates directed Weighted Geometric Chung-Lu model and divergence score for graph and embedding.

**Arguments**
* `edges::Array{Int,2}` array with edges definition (two whitespace separated vertices ids)
* `weights::Vector{Float64}` edges weights
* `comm::Array{Int,2}` assignment of vertices to communities
* `embed::Array{Float64,2}` array with vertices embeddings
* `distances::Vector{Float64}` distances between vertices
* `verbose::Bool` verbose switch, if true prints additional processing information
"""
function wGCL_directed(edges::Array{Int,2}, weights::Vector{Float64}, comm::Matrix{Int},
            embed::Matrix{Float64}, distances::Vector{Float64}, verbose::Bool=false)
    # default values
    delta = 0.001
    AlphaMax = 10.0
    AlphaStep = 0.25
    alpha_div_counter = 5
    auc_samples = 2000

    no_vertices = maximum(edges)
    no_edges = size(edges,1)
    verbose && println("$no_vertices vertices and $no_edges edges")

    # Read communities
    @assert size(comm,1) == no_vertices
    n_parts = maximum(comm)

    # Compute degrees
    degree_in = zeros(no_vertices)
    degree_out = zeros(no_vertices)
    star_check = zeros(Int64,no_vertices)
    for i in 1:no_edges
        a = weights[i]
        v1 = edges[i,1]
        v2 = edges[i,2]
        degree_out[v1]+=a
        degree_in[v2]+=a
        star_check[v1]+=1
        star_check[v2]+=1
    end

    # Check for edge case - star graph
    is_star = false
    # Star based on either in or out edges
    if !isnothing(findfirst(==(no_vertices-1),star_check)) && sum(star_check) == 2*(no_vertices-1)
        is_star = true
    # Star based on both in and out edges
    elseif !isnothing(findfirst(==(2*(no_vertices-1)),star_check)) && sum(star_check .== 2) == no_vertices-1
        is_star = true
    end
    verbose && is_star && println("Graph is a star in respect to either in or out edges")

    if is_star
        return [-1.0,zeros(5)...]
    end
    # Compute C-vector
    vect_len = Int(n_parts*n_parts)
    vect_C = zeros(Float64, vect_len)
    vect_B = zeros(Float64, vect_len)

    for i in 1:no_edges
        j = comm[edges[i,1]]
        k = comm[edges[i,2]]
        vect_C[(j-1)*n_parts + k] += weights[i]
    end

    # Indicator - internal to a community
    vect_I = falses(vect_len)
    for i in 1:(n_parts+1):vect_len
        vect_I[i] = true
    end
    best_div = best_div_ext = best_div_int = best_auc_err = best_auc = typemax(Float64)
    best_alpha = -1.0
    dim = size(embed,2)
    verbose && println("Embedding has $dim dimensions")
    # Loop over alpha's

    # Compute Euclidean distance vector D[] given embed and alpha
    p_len = no_vertices*(no_vertices+1) ÷ 2
    D = zeros(Float64, p_len)
    #Read distances
    @assert size(distances,1) == no_vertices
    for i in 1:no_vertices
        for j in i:no_vertices
            l = idx(no_vertices, i, j)
            if i == j
                D[l] = distances[i]
            else
                D[l] = dist(i,j,embed)
            end
        end
    end
    lo,hi = extrema(D)
    # Loop here - exclude Alpha = 0
    Tin = ones(no_vertices)
    Tin[degree_in.==0] .= 0.0
    Tout = ones(no_vertices)
    Tout[degree_out.==0] .= 0.0

    #Generate arrays of edges and non-edges
    NE = Tuple{Int64,Int64}[]
    for i in 1:no_vertices
        for j in 1:no_vertices
            if i!=j
                push!(NE,(i,j))
            end
        end
    end

    ## tuples of edges
    E = Tuple{Int64,Int64}[]
    for e in eachrow(edges)
        push!(E,tuple(e...))
    end

    ## tuples of non-edges
    NE = collect(setdiff(Set(NE),Set(E)));

    for alpha in AlphaStep:AlphaStep:(AlphaMax+delta)
        # Apply kernel (g(dist))
        GD = zeros(Float64, p_len)
        for i in 1:no_vertices
            for j in i:no_vertices
                k = idx(no_vertices,i,j)
                GD[k] = (D[k]-lo)/(hi-lo) # normalize to [0,1]
                GD[k] = (1-GD[k])^alpha # transform w.r.t. alpha
            end
        end
        # Learn GCL model numerically
        diff = 1.0
        epsilon = 0.9
        while diff > delta # stopping criterion
            Sin = zeros(no_vertices)
            Sout = zeros(no_vertices)
            for i in 1:no_vertices
                for j in i:no_vertices
                    k = idx(no_vertices, i, j)
                    tmp1 = Tin[i]*Tout[j]*GD[k]
                    tmp2 = Tin[j]*Tout[i]*GD[k]
                    Sin[i] += tmp1
                    Sin[j] += tmp2
                    Sout[i] += tmp2
                    Sout[j] += tmp1
                end
            end

            f = 0.0
            for i in 1:no_vertices
                if degree_in[i]>0
                    Tin[i] += epsilon*Tin[i]*(degree_in[i]/Sin[i]-1.0)
                    f = max(f,abs(degree_in[i]-Sin[i])) # convergence w.r.t. degrees
                end
                if degree_out[i]>0
                    Tout[i] += epsilon*Tout[i]*(degree_out[i]/Sout[i]-1.0)
                    f = max(f,abs(degree_out[i]-Sout[i])) # convergence w.r.t. degrees
                end
            end
            if f > diff
                epsilon *= 0.99
            end
            diff = f
            verbose && println("diff= $diff")
        end
        # Compute probas P[]
        P = zeros(Float64,2*p_len)
        for i in 1:no_vertices
            for j in 1:no_vertices
                a = min(i,j)
                b = max(i,j)
                P[no_vertices*(i-1)+j] = Tin[i]*Tout[j]*GD[idx(no_vertices,a,b)]
            end
        end
        ## random positive cases
        pos = [P[(no_vertices-1)*(E[i][1]-1)+E[i][2]-Int(E[i][2]>E[i][1])] for i in sample(1:length(E),auc_samples)]
        ## random negative cases
        neg = [P[(no_vertices-1)*(NE[i][1]-1)+NE[i][2]-Int(NE[i][2]>NE[i][1])] for i in sample(1:length(NE),auc_samples)]
        auc = 1 - sum(pos .> neg)/auc_samples
        ## auc estimate
        err = 1.96 * sqrt(auc*(1-auc)/auc_samples) ## error from 95% CI
        if auc < best_auc
            best_auc = auc
            best_auc_err = err
        end
        # Compute B-vector given P[] and comm[]
        vect_B = zeros(vect_len)
        for i in 1:no_vertices
            for j in 1:no_vertices
                m = (comm[i]-1)*n_parts + comm[j]
                vect_B[m] += P[no_vertices*(i-1)+j]
            end
        end
        x = JS(vect_C, vect_B, vect_I, true)
        y = JS(vect_C, vect_B, vect_I, false)
        f = (x+y)/2.0
        if f < best_div
            best_div = f
            best_alpha = alpha
            best_div_ext = x
            best_div_int = y
        else
            alpha_div_counter -= 1
            alpha_div_counter == 0 && break
        end
    end
    return [best_alpha, best_div, best_div_ext, best_div_int, best_auc, best_auc_err]
end
