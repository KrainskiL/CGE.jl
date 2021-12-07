var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"CurrentModule = CGE\nDocTestSetup = quote\n    using CGE\nend","category":"page"},{"location":"reference/#Landmarks","page":"Reference","title":"Landmarks","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"landmarks\nrunsplit\nsplit_cluster_diameter\nsplit_cluster_size\nsplit_cluster_rss\nsplit_cluster_rss2","category":"page"},{"location":"reference/#CGE.landmarks","page":"Reference","title":"CGE.landmarks","text":"landmarks(edges::Array{Int,2},vweights::Vector{Float64}, clusters::Vector{Vector{Int}},\n    embedding::Array{Float64,2}, verbose::Bool, land::Int, forced::Int, method::Function)\n\nArguments\n\nedges::Array{Int,2} array with edges definition (two whitespace separated vertices ids)\nweights::Vector{Float64} edges weights\nvweights::Vector{Float64} vertices weights\nclusters::Vector{Vector{Int}} vector of vectors indicating an initial 1-based cluster\n\nassignment of vertices.\n\nembedding::Array{Float64,2} array with vertices embeddings\ncomm::Array{Int,2} assignment of vertices to communities\nverbose::Bool verbose switch, if true prints additional processing information\nland::Int number of landmarks to generate\nforced::Int required maximum number of forced splits of a cluster\nmethod::Function method used to generate landmarks\ndirected::Bool flag for directed version of landmark-based graph creation\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.runsplit","page":"Reference","title":"CGE.runsplit","text":"runsplit(embedding, w, initial_clusters, n, s, rule)\n\nTake embedding with weights w where each column is a single observation and initial_clusters is a vector of vectors indicating an initial 1-based cluster assignment of vertices. Return a vector of 0-based assignments of vertices to n groups. Each cluster is guaranteed to be split to at most s landmarks (may be less if its size is less than s).\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.split_cluster_diameter","page":"Reference","title":"CGE.split_cluster_diameter","text":"split_cluster_diameter(m, w)\n\nSplits cluster m with weights w into two clusters along its first principal component so as to make both clusters have approximately the same diameter along the first principal component.\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.split_cluster_size","page":"Reference","title":"CGE.split_cluster_size","text":"split_cluster_size(m, w)\n\nSplits cluster m with weights w into two clusters along its first principal component so as to make the clusters have equal size (number of edges).\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.split_cluster_rss","page":"Reference","title":"CGE.split_cluster_rss","text":"split_cluster_rss(m, w)\n\nSplits cluster m with weights w into two clusters along its first principal component so as to minimize maximum RSS of one of the resulting clusters.\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.split_cluster_rss2","page":"Reference","title":"CGE.split_cluster_rss2","text":"split_cluster_rss2(m, w)\n\nSplits cluster m with weights w into two clusters along its first principal component so as to minimize maximum RSS of one of the resulting clusters.\n\nThis is a second slower version of the algorithm using sorting. Retained for testing purposes.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Divergence","page":"Reference","title":"Divergence","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"wGCL\nwGCL_directed","category":"page"},{"location":"reference/#CGE.wGCL","page":"Reference","title":"CGE.wGCL","text":"wGCL(edges::Array{Int,2}, weights::Vector{Float64}, comm::Array{Int,2},             embed::Array{Float64,2}, distances::Vector{Float64}, verbose::Bool = false)\n\nCalculates Weighted Geometric Chung-Lu model and divergence score for graph and embedding.\n\nArguments\n\nedges::Array{Int,2} array with edges definition (two whitespace separated vertices ids)\neweights::Vector{Float64} edges weights\ncomm::Array{Int,2} assignment of vertices to communities\nembed::Array{Float64,2} array with vertices embeddings\ndistances::Vector{Float64} distances between vertices\nvweights::Vector{Float64} landmarks total weights - used only with landmarks approximation\ninit_vweights::Vector{Float64} vector with original (full) vertices weights - used only with landmarks approximation\nv_to_l::Vector{Int} mapping from vertices to landmarks (landmarks membership) - used only with landmarks approximation\ninit_edges::Array{Int,2} array with original (full) graph edges - used only with landmarks approximation\ninit_eweights::Vector{Float64} vector with original (full) edges weights - used only with landmarks approximation\ninit_embed::Matrix{Float64} array with original embedding for full graph - used only with landmarks approximation\nsplit::Bool indicator for splitting JS divergence score (global score)\nseed::Int RNG seed for local measure score\nauc_samples::Int no. samples for local measure score\nverbose::Bool verbose switch, if true prints additional processing information\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.wGCL_directed","page":"Reference","title":"CGE.wGCL_directed","text":"wGCL_directed(edges::Array{Int,2}, weights::Vector{Float64}, comm::Array{Int,2},             embed::Array{Float64,2}, distances::Vector{Float64}, verbose::Bool = false)\n\nCalculates directed Weighted Geometric Chung-Lu model and divergence score for graph and embedding.\n\nArguments\n\nedges::Array{Int,2} array with edges definition (two whitespace separated vertices ids)\neweights::Vector{Float64} edges weights\ncomm::Array{Int,2} assignment of vertices to communities\nembed::Array{Float64,2} array with vertices embeddings\ndistances::Vector{Float64} distances between vertices\nvweights::Vector{Float64} landmarks total weights - used only with landmarks approximation\ninit_vweights::Vector{Float64} vector with original (full) vertices weights - used only with landmarks approximation\nv_to_l::Vector{Int} mapping from vertices to landmarks (landmarks membership) - used only with landmarks approximation\ninit_edges::Array{Int,2} array with original (full) graph edges - used only with landmarks approximation\ninit_eweights::Vector{Float64} vector with original (full) edges weights - used only with landmarks approximation\ninit_embed::Matrix{Float64} array with original embedding for full graph - used only with landmarks approximation\nsplit::Bool indicator for splitting JS divergence score (global score)\nseed::Int RNG seed for local measure score\nauc_samples::Int no. samples for local measure score\nverbose::Bool verbose switch, if true prints additional processing information\n\n\n\n\n\n","category":"function"},{"location":"reference/#Auxilary","page":"Reference","title":"Auxilary","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"dist\nJS\nidx","category":"page"},{"location":"reference/#CGE.dist","page":"Reference","title":"CGE.dist","text":"dist(i::Int, j::Int, embed::Array{Float64,2})\n\nCalculates Euclidian distance between two vectors from embedding array.\n\nArguments\n\nv1::Int index of first vector\nv2::Int index of second vector\nembed::Array{Float64,2} graph embedding array\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.JS","page":"Reference","title":"CGE.JS","text":"JS(vC::Vector{Float64}, vB::Vector{Float64},\n    vI::Vector{Int}, internal::Int, vLen::Int)\n\nJensen-Shannon divergence with Dirichlet-like prior.\n\nArguments\n\nvC::Vector{Float64} first distribution of edges within and between communities\nvB::Vector{Float64} second distribution of edges within and between communities\nvI::Vector{Int} indicator of internal (1) and external (0) edges w.r.t. communities, if empty compute overall JS distance\ninternal::Int internal JS distance switch, if 1 return internal, else return external\n\n\n\n\n\n","category":"function"},{"location":"reference/#CGE.idx","page":"Reference","title":"CGE.idx","text":"Change 2-dimensional into 1-dim index\n\n\n\n\n\n","category":"function"},{"location":"reference/#Clustering","page":"Reference","title":"Clustering","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"louvain_clust","category":"page"},{"location":"reference/#CGE.louvain_clust","page":"Reference","title":"CGE.louvain_clust","text":"louvain_clust(edges::String)\n\nCalculate communities in graph using Louvain algoritm\n\nArguments\n\nedges::String name of file with edges definition\n\n\n\n\n\nlouvain_clust(filename::String, edges::Array{Int,2}, weights::Array{Float64,1}))\n\nCalculate communities in weighted graph using Louvain algoritm\n\nArguments\n\nfilename::String name of file with edges definition\nedges::Array{Int,2} list of edges\nweights::Array{Float64,1} array of edges' weights\n\n\n\n\n\n","category":"function"},{"location":"#CGE.jl","page":"CGE.jl","title":"CGE.jl","text":"","category":"section"}]
}