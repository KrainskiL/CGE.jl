# CGE.jl
Julia package to compare graph embeddings.

| **Documentation** | **Build Status** |
|---------------|--------------|
|[![][docs-latest-img]][docs-dev-url]| [![Build Status][travis-img]][travis-url]  [![Coverage Status][codecov-img]][codecov-url] <br/>|

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-url]: https://KrainskiL.github.io/CGE.jl/dev
[docs-stable-url]: https://KrainskiL.github.io/CGE.jl/stable

[travis-img]: https://travis-ci.org/KrainskiL/CGE.jl.svg?branch=master
[travis-url]: https://travis-ci.org/KrainskiL/CGE.jl

[codecov-img]: https://codecov.io/gh/KrainskiL/CGE.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/KrainskiL/CGE.jl

## Details of the framework

Presented at [WAW2020](https://math.ryerson.ca/waw2020/) with publication in [Springer LNCS](https://www.springer.com/gp/book/9783030484774).

Detailed information can be found in the [paper](https://math.ryerson.ca/~pralat/papers/2020_WAW-Scalable_Embeddings.pdf).

Framework version without landmarks (written in C) is available under: https://github.com/ftheberge/Comparing_Graph_Embeddings

## Installation

The current version uses Julia 1.4. Install `CGE.jl` by running Julia, switching to package manager by pressing `]` and running:
```
add https://github.com/KrainskiL/CGE.jl
```

In order to use the CLI for the package you should locate CGE_CLI.jl file after the installation.

The directory with CGE_CLI.jl file and example files can be found by running the following command in Julia:
```julia
using CGE; cd(pwd, joinpath(dirname(pathof(CGE)), "..", "example"))
```
Make sure to copy the CLI file from this location (as it is read only).

Alternatively you can just download CGE_CLI.jl from GitHub repository. It is located in `example/` folder.

Finally you might also download the whole repository and extract the CGE_CLI.jl file from it.
```shell
git clone https://github.com/KrainskiL/CGE.jl
mv CGE.jl/example/CGE_CLI.jl .
julia CGE_CLI.jl
```

# Running the code

Code computes the Jenssen-Shannon divergence between two edge distributions:
(1) first one is based on the supplied graph clustering
(2) second one is based on the supplied embedding and a Geometric Chung-Lu model
When comparing embeddings, lower divergence is better.

Format:

```
julia CGE_CLI.jl -g edgelist_file -e embedding_file [-c clusters_file] [-v] [-l landmarks -f forced -m method]

## required flags:
-g: the edgelist (1 per line, whitespace separated, optionally with weights)
-e: the embedding (two formats accepted, see details below)
## optional flags:
-c: the communities (in vertices order, 1 per line), if not given calculated using Louvain algorithm
-v: verbose, printing additional information
-l: number of landmarks to create
-f: number of forced landmarks to be created
-m: chosen ladnmark creation method: `rss`, `rss2`, `size`, `diamater`
```

For instance, while in `example` folder run:

```julia
julia CGE_CLI.jl -g 100k.edgelist -c 100k.ecg -e 100k.embedding -l 200 -f 0 -m diameter
```
Result consists of 4 elements:
1. Best alpha
2. **Best divergence score**
3. Best divergence external score
4. Best divergence internal score
```
[0.25, 0.01483964683262605, 0.026810577776668364, 0.002868715888583737]
```
# File Formats

For a graph with `n` nodes, the nodes can be represented with numbers 1 to n or 0 to n-1.

Two input files are required to run the algorithm:
1. the undirected graph, represented by a sequence of edges, 1 per line and with optional weights in third column
2. the node embedding in on of the supported formats (see below)

Additional file with the node's cluster number (1 per line) may be provided. If it's missing communities are calculated automatically with Louvain algorithm.

## Example of graph (edgelist) file

Nodes can be 0-based or 1-based.
One edge per line with whitespace between nodes.

```
1 32
1 22
1 20
1 18
1 14
1 13
1 12
1 11
1 9
1 8
...
```

Additional weights may be provided in third column (both integers and floats are supported).

```
1 32 1.23
1 22 2.13
1 20 3.12
1 18 1.23
1 14 1.23
1 13 0.45
1 12 0.21
1 11 1.5
1 9 0.61
1 8 1.23
...
```

## Example of clustering file

Clusters can be 0-based or 1-based.
If not provided, clusters will be automatically calculated with Louvain algorithm.

First variant with clusters IDs only - must be ordered by nodes IDs

```
1
1
1
1
0
0
0
1
3
1
...
```

Second variant with clusters IDs and nodes IDs

```
1 1
2 1
3 1
4 1
5 0
6 0
8 0 
9 1
7 3
11 1
...
```
## Example of embedding file

Nodes are 0-based or 1-based in any order.
Three formats are supported

**First format - unordered embedding**

First column indicates node number, the rest of the line is d-dimensional embedding.

```
21 0.960689 -2.28209 3.65194 0.272646 -3.01281 1.0245 -0.329389 -2.95956
33 0.702187 -2.14331 4.25541 0.372346 -3.16427 1.41296 -0.390471 -4.49782
3 0.854487 -2.30527 4.10575 0.370613 -3.04878 1.46481 -0.120326 -4.02328
29 0.673825 -2.19518 4.00447 0.650003 -2.74663 0.757385 -0.505723 -3.2947
32 0.750248 -2.26306 4.04495 0.143616 -3.02735 1.49937 -0.400896 -4.04177
25 0.831608 -2.191 4.04712 0.786012 -2.85804 1.11308 -0.391722 -3.4645
28 1.14632 -2.20708 4.11004 0.338067 -2.86409 1.01202 -0.485711 -3.50161
...
```

**Second format - ordered embedding**

Only d-dimensional embedding in order of nodes is stored in file.

```
0.854487 -2.30527 4.10575 0.370613 -3.04878 1.46481 -0.120326 -4.02328
0.960689 -2.28209 3.65194 0.272646 -3.01281 1.0245 -0.329389 -2.95956
0.831608 -2.191 4.04712 0.786012 -2.85804 1.11308 -0.391722 -3.4645
1.14632 -2.20708 4.11004 0.338067 -2.86409 1.01202 -0.485711 -3.50161
0.702187 -2.14331 4.25541 0.372346 -3.16427 1.41296 -0.390471 -4.49782
0.673825 -2.19518 4.00447 0.650003 -2.74663 0.757385 -0.505723 -3.2947
0.750248 -2.26306 4.04495 0.143616 -3.02735 1.49937 -0.400896 -4.04177
...
```

**Third format - node2vec format**

First line contains number of nodes and dimension of the embedding. It's stripped during parsing and the rest of the file is handled as either first or second format.
```
500 8
21 0.960689 -2.28209 3.65194 0.272646 -3.01281 1.0245 -0.329389 -2.95956
33 0.702187 -2.14331 4.25541 0.372346 -3.16427 1.41296 -0.390471 -4.49782
3 0.854487 -2.30527 4.10575 0.370613 -3.04878 1.46481 -0.120326 -4.02328
29 0.673825 -2.19518 4.00447 0.650003 -2.74663 0.757385 -0.505723 -3.2947
32 0.750248 -2.26306 4.04495 0.143616 -3.02735 1.49937 -0.400896 -4.04177
25 0.831608 -2.191 4.04712 0.786012 -2.85804 1.11308 -0.391722 -3.4645
28 1.14632 -2.20708 4.11004 0.338067 -2.86409 1.01202 -0.485711 -3.50161
...
```