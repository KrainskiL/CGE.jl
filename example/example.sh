echo "Approximate algorithm with 200 landmarks"
julia ./CGE_CLI.jl -g 10k.edgelist -c 10k.ecg -e 10k.embedding -l 200 --seed 42
echo "Approximate algorithm with default number of landmarks"
julia ./CGE_CLI.jl -g 10k.edgelist -c 10k.ecg -e 10k.embedding -l --seed 42
echo "Approximate algorithm with 5 force splits of communities into landmarks"
julia ./CGE_CLI.jl -g 10k.edgelist -c 10k.ecg -e 10k.embedding -f 5 --seed 42
