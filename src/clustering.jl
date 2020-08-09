#######################
# Community detection #
#######################

"""
    louvain(edges::String)

Calculate communities in graph using [Louvain algoritm](https://sites.google.com/site/findcommunities/)

**Arguments**
* `edges::String` name of file with edges definition
"""
function louvain(edges::String)
    BASEPATH = dirname(dirname(pathof(CGE)))
    TMP = Sys.iswindows() ? ENV["TEMP"] : "/tmp"
    if Sys.iswindows()
        TMP = ENV["TEMP"]
        run(`$BASEPATH/bin/win/convert.exe -i $edges -o $TMP/louvain.bin`)
        run(pipeline(`$BASEPATH/bin/win/louvain.exe $TMP/louvain.bin -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/win/hierarchy.exe $TMP/louvain.txt -l 1`,"$edges.ecg"))
    else
        chmod(joinpath(BASEPATH,"bin/unix/"),0o744,recursive=true)
        run(`$BASEPATH/bin/unix/convert -i $edges -o $TMP/louvain.bin`)
        run(pipeline(`$BASEPATH/bin/unix/louvain $TMP/louvain.bin -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/unix/hierarchy $TMP/louvain.txt -l 1`,"$edges.ecg"))
    end
end

"""
    louvain(filename::String, edges::Array{Int,2}, weights::Array{Float64,1}))

Calculate communities in weighted graph using [Louvain algoritm](https://sites.google.com/site/findcommunities/)

**Arguments**
* `filename::String` name of file with edges definition
* `edges::Array{Int,2}` list of edges
* `weights::Array{Float64,1}` array of egdges' weights
"""
function louvain(filename::String, edges::Array{Int,2}, weights::Array{Float64,1})
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
    if Sys.iswindows()
        run(`$BASEPATH/bin/win/convert.exe -i $tmp_edges -o $TMP/louvain.bin -w $tmp_weights`)
        run(pipeline(`$BASEPATH/bin/win/louvain.exe $TMP/louvain.bin -w $tmp_weights -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/win/hierarchy.exe $TMP/louvain.txt -l 1`,"$filename.ecg"))
    else
        chmod(joinpath(BASEPATH,"bin/unix/"),0o744,recursive=true)
        run(`$BASEPATH/bin/unix/convert -i $tmp_edges -o $TMP/louvain.bin -w $tmp_weights`)
        run(pipeline(`$BASEPATH/bin/unix/louvain $TMP/louvain.bin -w $tmp_weights -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/unix/hierarchy $TMP/louvain.txt -l 1`,"$filename.ecg"))
    end
end