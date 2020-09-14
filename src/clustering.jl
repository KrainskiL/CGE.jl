#######################
# Community detection #
#######################

"""
louvain_clust(edges::String)

Calculate communities in graph using [Louvain algoritm](https://sites.google.com/site/findcommunities/)

**Arguments**
* `edges::String` name of file with edges definition
"""
function louvain_clust(edges::String)
    BASEPATH = dirname(dirname(pathof(CGE)))
    TMP = Sys.iswindows() ? ENV["TEMP"] : "/tmp"
    convertlouvain() do exe
        run(`$exe -i $edges -o $TMP/louvain.bin`)
    end
    louvain() do exe
        run(pipeline(`$exe $TMP/louvain.bin -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
    end
    hierarchy() do exe
        run(pipeline(`$exe $TMP/louvain.txt -l 1`,"$edges.ecg"))
    end
end

"""
louvain_clust(filename::String, edges::Array{Int,2}, weights::Array{Float64,1}))

Calculate communities in weighted graph using [Louvain algoritm](https://sites.google.com/site/findcommunities/)

**Arguments**
* `filename::String` name of file with edges definition
* `edges::Array{Int,2}` list of edges
* `weights::Array{Float64,1}` array of edges' weights
"""
function louvain_clust(filename::String, edges::Array{Int,2}, weights::Array{Float64,1})
    BASEPATH = dirname(dirname(pathof(CGE)))
    TMP = Sys.iswindows() ? ENV["TEMP"] : "/tmp"
    tmp_edges = joinpath(TMP,filename)
    open(tmp_edges, "w") do f
        for i in edges println(f, i) end
    end
    tmp_weights = joinpath(TMP,filename*".weights")
    open(tmp_weights, "w") do f
        for i in weights println(f, i) end
    end
    convertlouvain() do exe
        run(`$exe -i $tmp_edges -o $TMP/louvain.bin -w $tmp_weights`)
    end
    louvain() do exe
        run(pipeline(`$exe $TMP/louvain.bin -w $tmp_weights -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
    end
    hierarchy() do exe
        run(pipeline(`$exe $TMP/louvain.txt -l 1`,"$filename.ecg"))
    end
end
