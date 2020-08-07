#######################
# Community detection #
#######################

function louvain(edges::String)
    BASEPATH = dirname(dirname(pathof(CGE)))
    run(`$BASEPATH/bin/convert -i $edges -o /tmp/louvain.bin`)
    run(pipeline(`$BASEPATH/bin/louvain /tmp/louvain.bin -l -1 -q id_qual`, stdout="/tmp/louvain.txt"))
    run(pipeline(`$BASEPATH/bin/hierarchy /tmp/louvain.txt -l 1`,`awk '{print $2}'`,`tail -n +2`,"$edges.ecg"))
end