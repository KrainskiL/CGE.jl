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
        l = idx(n_parts, j, k)
        vect_C[l] += a
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
    for i in 1:no_vertices
        for j in i:no_vertices
            f = dist(i,j,embed)
            k = idx(no_vertices, i, j)
            D[k] = f
        end
    end
    lo1, hi1 = extrema(D)
    # Read distances
    @assert size(distances,1) == no_vertices
    for i in 1:no_vertices
        f = distances[i]
        k = idx(no_vertices, i, i)
        D[k] = f
    end
    lo2, hi2 = extrema(D)
    lo,hi = extrema([lo1,lo2,hi1,hi2])

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
            k = 0
            verbose && println("diff1 = $diff")
            for i in 1:no_vertices
                for j in i:no_vertices
                    k = idx(no_vertices, i, j)
                    tmp = T[i]*T[j]*GD[k]
                    S[i] += tmp
                    if i!=j S[j]+=tmp end
                end
            end
            verbose && println("diff2 = $diff")
            f = 0.0
            for i in 1:no_vertices
                move = epsilon*T[i]*(degree[i]/S[i]-1.0)
                T[i]+=move
                f = max(f,abs(degree[i]-S[i])) # convergence w.r.t. degrees
            end
            diff = f
        end
        # Compute probas P[]
        P = zeros(Float64,p_len)
        for i in 1:no_vertices
            for j in i:no_vertices
                k = idx(no_vertices,i,j)
                P[k] = T[i]*T[j]*GD[k]
            end
        end

        # Compute B-vector given P[] and comm[]
        vect_B = zeros(vect_len)
        for i in 1:no_vertices
            for j in i:no_vertices
                k,l = extrema([comm[i],comm[j]])
                m = Int(n_parts*(k-1)-(k-1)*(k-2)/2+l-k+1)
                vect_B[m] += P[idx(no_vertices,i,j)]
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
    epsilon = 0.1
    delta = 0.001
    AlphaMax = 10.0
    AlphaStep = 0.25

    no_vertices = maximum(edges)
    no_edges = size(edges,1)
    verbose && println("$no_vertices vertices and $no_edges edges")

    # Read communities
    @assert size(comm,1) == no_vertices
    n_parts = maximum(comm)

    # Compute degrees
    degree_in = zeros(no_vertices)
    degree_out = zeros(no_vertices)
    for i in 1:no_edges
        a = weights[i]
        degree_out[edges[i,1]]+=a
        degree_in[edges[i,2]]+=a
    end

    # Compute C-vector
    vect_len = Int(n_parts*n_parts)
    vect_C = zeros(Float64, vect_len)
    vect_B = zeros(Float64, vect_len)
 
    for i in 1:no_edges
        a = weights[i]
        j = comm[edges[i,1]]
        k = comm[edges[i,2]]
        l = (j-1)*n_parts + k
        vect_C[l] += a
    end

    # Indicator - internal to a community
    vect_I = falses(vect_len)
    j=1
    for i in 1:n_parts
        vect_I[j] = true
        j+=(n_parts+1)
    end
    best_div = best_div_ext = best_div_int = typemax(Float64)
    best_alpha = -1.0
    dim = size(embed,2)
    verbose && println("Embedding has $dim dimensions")

    # Loop over alpha's

    # Compute Euclidean distance vector D[] given embed and alpha
    p_len = no_vertices*(no_vertices+1) ÷ 2
    D = zeros(Float64, p_len)
    for i in 1:no_vertices
        for j in i:no_vertices
            f = dist(i,j,embed)
            k = idx(no_vertices, i, j)
            D[k] = f
        end
    end
    lo1, hi1 = extrema(D)
    #Read distances
    @assert size(distances,1) == no_vertices
    for i in 1:no_vertices
        f = distances[i]
        k = idx(no_vertices, i, i)
        D[k] = f
    end
    lo2, hi2 = extrema(D)
    lo,hi = extrema([lo1,lo2,hi1,hi2])
    # lo,hi = extrema(D)
    # Loop here - exclude Alpha = 0
 
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
        Tin = ones(no_vertices)
        Tin[degree_in.==0] .= 0.0
        Tout = ones(no_vertices)
        Tout[degree_in.==0] .= 0.0   
        diff = 1.0

        while diff > delta # stopping criterion
            Sin = zeros(no_vertices)
            Sout = zeros(no_vertices)
            verbose && println("diff1 = $diff")
            k=0
            for i in 1:no_vertices
                for j in i:no_vertices   
                    k = idx(no_vertices, i, j)
                    Sin[i] += Tin[i]*Tout[j]*GD[k]
                    Sin[j] += Tin[j]*Tout[i]*GD[k]
                    Sout[i] += Tout[i]*Tin[j]*GD[k]
                    Sout[j] += Tout[j]*Tin[i]*GD[k]
                end
            end          
            
            f = 0.0
            for i in 1:no_vertices
                if degree_in[i]>0
                    move = epsilon*Tin[i]*(degree_in[i]/Sin[i]-1.0)
                    Tin[i]+=move
                    f = max(f,abs(degree_in[i]-Sin[i])) # convergence w.r.t. degrees
                end
                if degree_out[i]>0
                    move = epsilon*Tout[i]*(degree_out[i]/Sout[i]-1.0)
                    Tout[i]+=move
                    f = max(f,abs(degree_out[i]-Sout[i])) # convergence w.r.t. degrees
                end
            end
            diff = f
            verbose && println("diff2 = $diff")
        end

        # Compute probas P[]
        P = zeros(Float64,2*p_len)
        for i in 1:no_vertices
            for j in 1:no_vertices
                if i != j
                    k = (no_vertices-1)*(i-1)+j-Int(j>i)
                    a = min(i,j)
                    b = max(i,j)
                    l = idx(no_vertices,a,b)
                    P[k] = Tin[i]*Tout[j]*GD[l]
                end
            end
        end

        # Compute B-vector given P[] and comm[]
        vect_B = zeros(vect_len)
        for i in 1:no_vertices
            for j in 1:no_vertices
                if i != j
                    m = (comm[i]-1)*n_parts + comm[j]
                    vect_B[m] += P[(no_vertices-1)*(i-1)+j-Int(j>i)]
                end
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
        end
    end
    return [best_alpha, best_div, best_div_ext, best_div_int]
end