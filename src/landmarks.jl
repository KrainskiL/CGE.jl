##################
# Landmarks code #
##################

struct Landmark
    what::Vector{Int}
    value::Float64
end

# priority queue implementation to avoid DataStructures.jl dependency

function landmark_put!(pq::Vector{Landmark}, what::Vector{Int}, value::Float64)
    p = Landmark(what, value)
    push!(pq, p)
    i = length(pq)
    @inbounds while (j = div(i, 2)) >= 1
        if value < pq[j].value
            pq[i] = pq[j]
            i = j
        else
            break
        end
    end
    pq[i] = p
end

function landmark_pop!(pq::Vector{Landmark})
    x = pq[1]
    y = pop!(pq)
    if !isempty(pq)
        i = 1
        len = length(pq)
        @inbounds while (l = 2i) <= len
            r = 2i + 1
            j = r > len || pq[l].value < pq[r].value ? l : r
            if pq[j].value < y.value
                pq[i] = pq[j]
                i = j
            else
                break
            end
        end
        pq[i] = y
    end
    return x.what, x.value
end

# On-line SSE calculation, to OnlineStats.jl dependency

struct WSSE
    ss::Float64 # weighted sum of squared values
    s::Float64  # weighted sum of values
    ws::Float64 # sum of weights
end

WSSE(x::AbstractVector{<:Real}, w::AbstractVector{<:Real}) =
    if isempty(x)
        WSSE(0.0, 0.0, 0.0)
    else
        WSSE(sum(v->v[2]*v[1]^2, zip(x, w)), dot(w, x), sum(w))
    end
WSSE(x::Real, w::Real) = WSSE(w*x^2, w*x, w)

# this is numerically unstable but we are OK even if we have some inaccuracy
wsse(x::WSSE) = x.ss - x.s^2 / x.ws
Base.:+(x1::WSSE, x2::WSSE) = WSSE(x1.ss + x2.ss, x1.s + x2.s, x1.ws + x2.ws)
Base.:-(x1::WSSE, x2::WSSE) = WSSE(x1.ss - x2.ss, x1.s - x2.s, x1.ws - x2.ws)

# weighted mean to avoid StatsBase.jl dependency

function matrix_w_mean(m, w)
    sw = sum(w)
    res = zeros(1, size(m, 2))
    for i in axes(m, 2)
        for j in axes(m, 1)
            @inbounds res[1, i] += m[j, i] * w[j]
        end
        res[1, i] /= sw
    end
    return res
end

"""
    split_cluster_rss2(m, w)

Splits cluster `m` with weights `w` into two clusters along its first principal
component so as to minimize maximum RSS of one of the resulting clusters.

This is a second slower version of the algorithm using sorting.
Retained for testing purposes.
"""
function split_cluster_rss2(m, w)
    @assert size(m, 1) > 1
    if size(m, 1) == 2
        return [1], [2]
    end
    y = (m .- matrix_w_mean(m, w)) .* sqrt.(w)
    yᵀwy = transpose(y) * y
    z = y * (@view eigvecs(yᵀwy)[:, end])
    p = sortperm(z)
    low = 1
    high = length(p)
    rss_low = [WSSE(m[p[1], i], w[p[1]]) for i in axes(m, 2)]
    rss_high = [WSSE(m[p[end], i], w[p[end]]) for i in axes(m, 2)]
    while low + 1 < high
        if sum(wsse, rss_low) < sum(wsse, rss_high)
            low += 1
            for (i, v) in enumerate(rss_low)
                rss_low[i] = v + WSSE(m[p[low], i], w[p[low]])
            end
        else
            high -= 1
            for (i, v) in enumerate(rss_high)
                rss_high[i] = v + WSSE(m[p[high], i], w[p[high]])
            end
        end
    end
    moved_low = false
    while low > 1
        rss_low_tmp = [rss_low[i] - WSSE(m[p[low], i], w[p[low]]) for i in axes(rss_low, 1)]
        rss_high_tmp = [rss_high[i] + WSSE(m[p[low], i], w[p[low]]) for i in axes(rss_high, 1)]
        if max(sum(wsse, rss_low_tmp), sum(wsse, rss_high_tmp)) < max(sum(wsse, rss_low), sum(wsse, rss_high))
            moved_low = true
            low -= 1
            high -= 1
            rss_low = rss_low_tmp
            rss_high = rss_high_tmp
        else
            break
        end
    end
    if !moved_low
        while high < length(p)
            rss_low_tmp = [rss_low[i] + WSSE(m[p[high], i], w[p[high]]) for i in axes(rss_low, 1)]
            rss_high_tmp = [rss_high[i] - WSSE(m[p[high], i], w[p[high]]) for i in axes(rss_high, 1)]
            if max(sum(wsse, rss_low_tmp), sum(wsse, rss_high_tmp)) < max(sum(wsse, rss_low), sum(wsse, rss_high))
                low += 1
                high += 1
                rss_low = rss_low_tmp
                rss_high = rss_high_tmp
            else
                break
            end
        end
    end
    return view(p, 1:low), view(p, high:length(p))
end

"""
    split_cluster_rss(m, w)

Splits cluster `m` with weights `w` into two clusters along its first principal
component so as to minimize maximum RSS of one of the resulting clusters.
"""
function split_cluster_rss(m, w)
    @assert size(m, 1) > 1
    if size(m, 1) == 2
        return [1], [2]
    end
    y = (m .- matrix_w_mean(m, w)) .* sqrt.(w)
    yᵀwy = transpose(y) * y
    z = y * (@view eigvecs(yᵀwy)[:, end])
    l1 = Int[argmin(z)]
    l2 = Int[argmax(z)]
    if l1 == l2
        throw(ErrorException("Trying to split homogenous cluster"))
    end
    gray = setdiff(axes(z, 1), l1, l2)
    rss_low = [WSSE(m[l1[1],i]^2*w[l1[1]], m[l1[1],i]*w[l1[1]], w[l1[1]]) for i in axes(m, 2)]
    rss_high = [WSSE(m[l2[1],i]^2*w[l2[1]], m[l2[1],i]*w[l2[1]], w[l2[1]]) for i in axes(m, 2)]
    med = median(z)
    t1 = Int[]
    t2 = Int[]
    while true
        empty!(t1)
        empty!(t2)
        for i in gray
            if z[i] < med
                push!(t1, i)
            else
                push!(t2, i)
            end
        end
        rss_low_tmp = [wssei + WSSE(view(m, t1, i), view(w, t1)) for (i, wssei) in enumerate(rss_low)]
        rss_high_tmp = [wssei + WSSE(view(m, t2, i), view(w, t2)) for (i, wssei) in enumerate(rss_high)]
        if sum(wsse, rss_low_tmp) < sum(wsse, rss_high_tmp)
            length(t1) == 0 && break
            rss_low = rss_low_tmp
            append!(l1, t1)
            gray = copy(t2)
        else
            length(t2) == 0 && break
            rss_high = rss_high_tmp
            append!(l2, t2)
            gray = copy(t1)
        end
        isempty(gray) && break
        med = median(view(z, gray))
    end
    if length(gray) > 0
        rss_low_tmp = [wssei + WSSE(view(m, gray, i), view(w, gray)) for (i, wssei) in enumerate(rss_low)]
        rss_high_tmp = [wssei + WSSE(view(m, gray, i), view(w, gray)) for (i, wssei) in enumerate(rss_high)]
        if max(sum(wsse, rss_low_tmp), sum(wsse, rss_high)) < max(sum(wsse, rss_low), sum(wsse, rss_high_tmp))
            append!(l1, gray)
        else
            append!(l2, gray)
        end
    end
    return l1, l2
end

"""
    split_cluster_size(m, w)

Splits cluster `m` with weights `w` into two clusters along its first principal
component so as to make the clusters have equal size (number of edges).
"""
function split_cluster_size(m, w)
    @assert size(m, 1) > 1
    if size(m, 1) == 2
        return [1], [2]
    end
    y = (m .- matrix_w_mean(m, w)) .* sqrt.(w)
    yᵀwy = transpose(y) * y
    z = y * (@view eigvecs(yᵀwy)[:, end])
    med = median(z)

    low = Int[]
    high = Int[]
    for (i, c) in enumerate(z)
        if c == med
            push!(length(low) < length(high) ? low : high, i)
        else
            push!(c < med ? low : high, i)
        end
    end
    return low, high
end

"""
    split_cluster_diameter(m, w)

Splits cluster `m` with weights `w` into two clusters along its first principal
component so as to make both clusters have approximately the same diameter along
the first principal component.
"""
function split_cluster_diameter(m, w)
    @assert size(m, 2) > 1
    if size(m, 1) == 2
        return [1], [2]
    end
    y = (m .- matrix_w_mean(m, w)) .* sqrt.(w)
    yᵀwy = transpose(y) * y
    z = y * (@view eigvecs(yᵀwy)[:, end])
    avg = mean(extrema(z))

    low = Int[]
    high = Int[]
    for (i, c) in enumerate(z)
        if c == avg
            push!(length(low) < length(high) ? low : high, i)
        else
            push!(c < avg ? low : high, i)
        end
    end
    return low, high
end

total_rss(m::AbstractMatrix, w) = sum(x -> wsse(WSSE(x, w)), eachcol(m))

"""
    runsplit(embedding, w, initial_clusters, n, s, rule)

Take `embedding` with weights `w` where each column is a single observation and
`initial_clusters` is a vector of vectors indicating an initial 1-based cluster
assignment of vertices. Return a vector of 0-based assignments of vertices to `n` groups.
Each cluster is guaranteed to be split to at most `s` landmarks (may be less if its size is less than `s`).
"""
function runsplit(embedding, w, initial_clusters, n, s, rule)
    landmarks = Vector{Landmark}()
    for cluster in sort(initial_clusters)
        if length(cluster) <= s
            for j in cluster
                landmark_put!(landmarks, [j], eps())
            end
        else
            locallandmarks = Vector{Landmark}()
            landmark_put!(locallandmarks, cluster, -total_rss(view(embedding, cluster, :), view(w, cluster)))
            while length(locallandmarks) < s
                idxs, sse = landmark_pop!(locallandmarks)
                low, high = rule(view(embedding, idxs, :), view(w, idxs))
                idxs_low = idxs[low]
                if length(idxs_low) > 1
                    landmark_put!(locallandmarks, idxs_low, -total_rss(view(embedding, idxs_low, :), view(w, idxs_low)))
                elseif length(idxs_low) == 1
                    landmark_put!(locallandmarks, idxs_low, eps()) # make sure 1-length cluster is in tail of the queue
                else
                    throw(ErrorException("Unexpected empty cluster generated"))
                end
                idxs_high = idxs[high]
                if length(idxs_high) > 1
                    landmark_put!(locallandmarks, idxs_high, -total_rss(view(embedding, idxs_high, :), view(w, idxs_high)))
                elseif length(idxs_high) == 1
                    landmark_put!(locallandmarks, idxs_high, eps()) # make sure 1-length cluster is in tail of the queue
                else
                    throw(ErrorException("Unexpected empty cluster generated"))
                end
            end
            while !isempty(locallandmarks)
                idxs, sse = landmark_pop!(locallandmarks)
                landmark_put!(landmarks, idxs, sse)
            end
        end
    end

    while length(landmarks) < n
        idxs, sse = landmark_pop!(landmarks)
        low, high = rule(view(embedding, idxs, :), view(w, idxs))
        idxs_low = idxs[low]
        if length(idxs_low) > 1
            landmark_put!(landmarks, idxs_low, -total_rss(view(embedding, idxs_low, :), view(w, idxs_low)))
        elseif length(idxs_low) == 1
            landmark_put!(landmarks, idxs_low, eps()) # make sure 1-length cluster is in tail of the queue
        else
            throw(ErrorException("Unexpected empty cluster generated"))
        end
        idxs_high = idxs[high]
        if length(idxs_high) > 1
            landmark_put!(landmarks, idxs_high, -total_rss(view(embedding, idxs_high, :), view(w, idxs_high)))
        elseif length(idxs_high) == 1
            landmark_put!(landmarks, idxs_high, eps()) # make sure 1-length cluster is in tail of the queue
        else
            throw(ErrorException("Unexpected empty cluster generated"))
        end
    end

    group_ids = fill(-1, size(embedding, 1))
    i = 0
    for vs in landmarks
        group_ids[vs.what] .= i
        i += 1
    end
    @assert all(>=(0), group_ids)
    return group_ids
end

"""
    landmarks(edges::Array{Int,2},vweights::Vector{Float64}, clusters::Vector{Vector{Int}},
        embedding::Array{Float64,2}, verbose::Bool, land::Int, forced::Int, method::Function)

**Arguments**
* `edges::Array{Int,2}` array with edges definition (two whitespace separated vertices ids)
* `vweights::Vector{Float64}` vertices degrees
* `clusters::Vector{Vector{Int}}` vector of vectors indicating an initial 1-based cluster
assignment of vertices.
* `embedding::Array{Float64,2}` array with vertices embeddings
* `comm::Array{Int,2}` assignment of vertices to communities
* `verbose::Bool` verbose switch, if true prints additional processing information
* `land::Int` number of landmarks to generate
* `forced::Int` required maximum number of forced splits of a cluster
* `method::Function` method used to generate landmarks
"""
function landmarks(edges::Array{Int,2}, weights::Vector{Float64}, vweights::Vector{Float64},
                    clusters::Vector{Vector{Int}}, comm::Array{Int,2},embedding::Array{Float64,2},
                    verbose::Bool, land::Int, forced::Int, method::Function)

    verbose && println("Starts landmark generation")
    rows_embed, dim = size(embedding)
    unique_rows = size(unique(embedding, dims=1), 1) #performance 0.5 sec for 100k vertices

    if land > unique_rows
        @warn "Requested number of clusters larger than unique no. embeddings. Truncating to $unique_rows landmarks."
        land = unique_rows
    end

    landmarks = runsplit(embedding, vweights, clusters, land, forced, method)
    landmarks .+=1

    verbose && println("Landmarks generated")

    # Count no. landmarks
    N = Int(maximum(landmarks))
    verbose && println("Using $N landmarks")

    embed = zeros(N,dim)
    lweight = zeros(N)

    # read from 0-based embedding file
    for i in 1:rows_embed
        what = landmarks[i]
        lweight[what]+=1
        for j in 1:dim
            embed[what,j]+=embedding[i,j]
        end
    end

    # normalize landmark based embedding
    for i in 1:N
        for j in 1:dim
            embed[i,j] /= lweight[i]
        end
    end

    # compute landmark internal distances d_ii
    dii = zeros(N)
    for i in 1:rows_embed
        what = landmarks[i]
        dist = 0
        for j in 1:dim
            dist += (embed[what,j]-embedding[i,j])^2
        end
        dii[what] += dist
    end

    # normalize
    for i in 1:N
        if lweight[i]>0
            dii[i] /= lweight[i]
            dii[i] = sqrt(dii[i])
        end
    end

    # re-write clusters for the landmarks
    cluster = zeros(Int,N)
    for i in 1:rows_embed
        cluster[landmarks[i]]=comm[i]
    end
    cluster = reshape(cluster,:,1)

    # Output weighted landmark_edgelist
    wedges = zeros(Float64,N,N)
    for i in 1:size(edges,1)
        a,b = extrema([landmarks[edges[i,1]],landmarks[edges[i,2]]])
        wedges[a,b] += weights[i]
    end

    # re-write
    landmark_edges = Array{Float64}(undef, Int(N*(N+1)/2), 3)
    for i in 1:N
        for j in i:N
            landmark_edges[idx(N,i,j),:] = [i j wedges[i,j]]
        end
    end
    landmark_edges = landmark_edges[landmark_edges[:,3] .>0,:]
    weights = landmark_edges[:,3]
    landmark_edges = Int.(landmark_edges[:,1:2])

    return dii, embed, cluster, landmark_edges, weights
end