#######################
# Community detection #
#######################

function louvain(edges::String)
    BASEPATH = dirname(dirname(pathof(CGE)))
    if Sys.iswindows()
        TMP = ENV["TEMP"]
        run(`$BASEPATH/bin/win/convert.exe -i $edges -o $TMP/louvain.bin`)
        run(pipeline(`$BASEPATH/bin/win/louvain.exe $TMP/louvain.bin -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/win/hierarchy.exe $TMP/louvain.txt -l 1`,"$edges.ecg"))
    else
        run(`$BASEPATH/bin/unix/convert -i $edges -o /tmp/louvain.bin`)
        run(pipeline(`$BASEPATH/bin/unix/louvain /tmp/louvain.bin -l -1 -q id_qual`, stdout="/tmp/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/unix/hierarchy /tmp/louvain.txt -l 1`,"$edges.ecg"))
    end
end

function louvain(filename::String, edges::Array{Int,2}, weights::Array{Float64,1})
    BASEPATH = dirname(dirname(pathof(CGE)))
    TMP = Sys.iswindows() ? ENV["TEMP"] : "/tmp"
    tmp_edges = TMP * filename 
    open(tmp_edges, "w") do f
        for i in edges println(f, i) end
    end
    tmp_weights = TMP * filename*".weights"
    open(tmp_weights, "w") do f
        for i in weights println(f, i) end
    end
    if Sys.iswindows()
        run(`$BASEPATH/bin/win/convert.exe -i $tmp_edges -o $TMP/louvain.bin -w $tmp_weights`)
        run(pipeline(`$BASEPATH/bin/win/louvain.exe $TMP/louvain.bin -w $tmp_weights -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/win/hierarchy.exe $TMP/louvain.txt -l 1`,"$filename.ecg"))
    else
        run(`$BASEPATH/bin/unix/convert -i $tmp_edges -o $TMP/louvain.bin -w $tmp_weights`)
        run(pipeline(`$BASEPATH/bin/unix/louvain $TMP/louvain.bin -w $tmp_weights -l -1 -q id_qual`, stdout="$TMP/louvain.txt"))
        run(pipeline(`$BASEPATH/bin/unix/hierarchy $TMP/louvain.txt -l 1`,"$filename.ecg"))
    end
end