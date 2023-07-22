# koala-networkit

This project is a KOALA library fork built on top of the structures provided by the NetworKit library.

#### Table of Contents
1. [Overview of the library](#overview)
    * [NetworKit](#networkit)
    * [KOALA](#koala)
    * [List of algorithms](#algorithms)
2. [Installation](#installation)
3. [Usage](#usage)
4. [References](#references)

## <a name="overview"></a>Overview of the library

<p align="justify">
In order to provide a simple, fast, powerful unified interface for a broad class of algorithms we joined the power and scalability of a fast, popular, and modern NetworKit library with an equally unique set of algorithmic tools for classical discrete optimization problems provided by the KOALA library.
</p>

### <a name="networkit"></a>NetworKit

<p align="justify">
<a href="https://networkit.github.io/">NetworKit</a> is an open-source general-purpose network analysis and graph mining library written in C++ (see <a href="#first">[1]</a> for details).
Graphs are represented in a very compact form while the efficiency of changes to their structure, such as node and edge additions or deletions, is preserved.
The memory-saving design is related to the main aim of the library: the analysis of large-scale random graphs and real-world networks.
</p>

<p align="justify">
The library contains algorithms mostly for the structural analysis of massive complex real-world networks, i.e. finding global network properties (such as density of a graph or its diameter), centrality measures (betweenness, closeness, local clustering), density and motifs measures (clustering coefficient, triangle counting, clique detection), and community detection. It also includes implementations of many generative network models such as Erdős–Rényi, Barabási–Albert or the hyperbolic unit-disk model.
It has some procedures related to the classical graph problems e.g. Edmonds-Karp maximum flow algorithm, but it is lacking depth in this dimension, for example, there are no algorithms for graph coloring.
</p>

<p align="justify">
The authors of this library emphasize the usage of parallelization, heuristics, and efficient data structures to deal with computationally intensive problems on large data.
For design, they aim at a modular architecture with encapsulation of algorithms into software components (classes and modules) for extensibility and code reuse. The code is written very clearly with promoting the compatibility of good coding practices and new standards of C++ in mind.
For users, the library provides seamless integration with Python and its libraries e.g. pandas, matplotlib. It also contains interfaces to some external products e.g. graph visualization tool Gephi.

<p align="justify">
This library comes out as one the fastest in the field (see <a href="https://www.timlrx.com/blog/benchmark-of-popular-graph-network-packages-v2">here</a> for a comparison), it is quite readable and efficient in practice, especially compared to ad-hoc made up solutions.
Although not as popular as <a href="https://github.com/snap-stanford/snap">Stanford Network Analysis Platform (SNAP)</a> library, and not as endowed by Chan Zuckerberg Initiative as <a href="https://github.com/igraph/igraph">igraph</a> library, NetworKit still has certain lively community (over 450 stars and 160 forks on GitHub) keeping with its development.
</p>

### <a name="koala"></a>Koala

<p align="justify">
<a href="http://web.archive.org/web/20200721235426/http://koala.os.niwa.gda.pl/api/description.html">KOALA</a> is an open-source library of C++ templates, developed at the Gdansk University of Technology, Department of Algorithms and System Modeling (see <a href="#second">[2]</a> for details).
Its main part consists of an implementation of a broad set of procedures in the fields of algorithmic graph theory and network problems in discrete optimization. In particular, the library contains algorithms for:
</p>

- graph coloring:
  - vertex coloring: heuristics (greedy, LF, SL, SLF, GIS, also with color interchange), $\Delta$-coloring, exact exponential-time,
  - edge coloring: heuristics (greedy, greedy with interchange), $(\mu + \Delta)$-coloring,
  - list vertex coloring, list edge coloring, interval vertex coloring, interval edge coloring,
- graph search:
  - graph traversal: DFS, BFS, lexicographic BFS (with pre- and post-order versions),
  - shortest paths: Dijkstra, Bellman-Ford, Johnson, Floyd-Warshall,
  - spanning forest: Kruskal, Prim,
  - computation of connected components and biconnected components, negative cycle detection,
- flow problems: maximum flow (Fulkerson-Ford, Dinic), minimum cost maximum flow, minimum edge or vertex cut, Gomory-Hu tree,
- independent set (heuristics, approximations, and exact exponential-time) and factorization,
- task scheduling on parallel identical processors with classical measures of the schedule cost: $C_{max}$, $L_{max}$, $\sum C_i$, $\sum U_i$, $\sum T_i$,
- graph recognition for graph families: empty graphs, cliques, paths, caterpillar, trees, forests, cycles, connected, complete $k$-partite, regular, subcubic, block, bipartite, chordal, comparability, interval, split, cographs,
- graph generation for several graph families (cliques, paths, cycles, fans, wheels, caterpillars, complete $k$-partite, regular trees) and random graphs models (Erdős–Rényi, Barabási–Albert, Watts–Strogatz),
- operation on graphs: closures (reflexive, symmetric, transitive), line graphs, products of graphs (Cartesian, tensor, lexicographic, strong).
- dominating set: exact exponential-time

<p align="justify">
The library was built on a custom, versatile, object-oriented templated graph structure, capable of handling edges of multiple types (undirected, directed, loops) at the same time.
The creators wanted to provide a convenient and friendly user experience through interface procedures on various levels of complexity.
At the most basic level, they presented generic versions of the algorithms, but using C++ language mechanisms they also allow more advanced programmers e.g. to customize the data structures (arrays, maps, priority queues) used by the algorithms instead of using the STL containers. For example, users can pick the structures which are tailored to guarantee computational complexity or the ones which are tailored to special practical cases.
</p>

<p align="justify">
Moreover, they set up an online graph editor <a href="https://stos.eti.pg.gda.pl/~kmocet/zgred/1.1.22/zgred.html">Zgred</a>, written in JavaScript. It allows to create, edit and visualize graphs. Furthermore, it is capable of running several algorithms from the library.
</p>

### <a name="algorithms"></a>List of algorithms

1. [Reading and writing graphs](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/io): [graph6](https://users.cecs.anu.edu.au/~bdm/data/formats.html), [sparse6](https://users.cecs.anu.edu.au/~bdm/data/formats.html), [digraph6](https://users.cecs.anu.edu.au/~bdm/data/formats.html),[DIMACS](http://prolland.free.fr/works/research/dsat/dimacs.html), [DIMACS binary](https://mat.tepper.cmu.edu/COLOR/format/README.binformat) formats
2. [Graph recognition](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/recognition/): [perfect graphs](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/recognition/PerfectGraphRecognition.hpp)
3. [Graph traversal](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/traversal/): [BFS](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/traversal/BFS.hpp), [DFS](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/traversal/DFS.hpp)
4. [Flow algorithms](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/flow/)
    1. [Maximum flow](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/flow/MaximumFlow.hpp): King-Rao-Tarjan
5. [Vertex coloring](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/coloring/)
    1. [Greedy heuristics](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/coloring/GreedyVertexColoring.hpp): RandomSequential, LargestFirst, SmallestLast, SaturatedLargestFirst, GreedyIndependentSet
    2. [Exact exponential-time algorithms](https://github.com/krzysztof-turowski/koala-networkit/blob/master/include/coloring/ExactVertexColoring.hpp): Brown, Christofides, Brélaz, Korman
    3. [Grötschel-Lovász-Schrijver algorithm for perfect graphs](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/coloring/PerfectGraphVertexColoring.hpp)
6. [Minimum dominating set](https://github.com/krzysztof-turowski/koala-networkit/tree/master/include/dominating_set/): Grandoni, Fomin-Grandoni-Kratsch, van Rooij-Bodlaender, Fomin-Kratsch-Woeginger, Schiermeyer
7. Minimum set cover: [exact branch and reduce](https://github.com/krzysztof-turowski/koala-networkit/blob/master/include/set_cover/BranchAndReduceMSC.hpp)


For further planned changes, see the [Issues](https://github.com/krzysztof-turowski/koala-networkit/issues/) section.

### <a name="datasets"></a>List of datasets

<p align="justify">
To assess the speed of the algorithms we use primarily the publicly available Stanford Large Network Dataset, a set of real-world networks (social, information, biological, etc.) and datasets <a href="#third">[3]</a>.
</p>

## <a name="installation"></a>Installation

```bash
cmake -B build
cmake --build build --parallel 4
```
> Note: it may take a while to download and compile dependencies (e.g. googletest, networkit, and boost).

## <a name="usage"></a>Usage

## <a name="references"></a>References

[[1]](https://www.cambridge.org/core/journals/network-science/article/networkit-a-tool-suite-for-largescale-complex-network-analysis/03DB673D73EDC84C0A143864FFA17831)<a name="first"></a> Christian Staudt, Aleksejs Sazonovs, Henning Meyerhenke, <i>NetworKit: A tool suite for large-scale complex network analysis</i>, Network Science 4(4):508-530, 2016.

[[2]](https://task.gda.pl/files/quart/TQ2015/04/tq419r-c.pdf)<a name="second"></a> Krzysztof Giaro, Krzysztof Ocetkiewicz, Tomasz Goluch, <i>KOALA Graph Theory Internet Service</i>, TASK Quarterly 19(4):455-470, 2015.

[[3]](http://snap.stanford.edu/data)<a name="third"></a> Jure Leskovec, Rok Sosic, <i>SNAP Datasets: Stanford Large Network Dataset Collection</i>, 2014.
