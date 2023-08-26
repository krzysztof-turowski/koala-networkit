/*
 * IndependentSet.cpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#include <independent_set/IndependentSet.hpp>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <limits>
#include <set>
#include <bitset>

namespace Koala {

void printGraph(const NetworKit::Graph* graph) { // TODO: its just debug
    graph->forNodes([&](NetworKit::node v) {
        std::cout << v << " ";
    });
    std::cout << std::endl;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        std::cout << "{" << u << "," << v << "}, ";
    });
    std::cout << std::endl;
}

std::vector<NetworKit::node> IndependentSet::getNeighbors(NetworKit::node v) const {
    std::vector<NetworKit::node> neighbors;
    graph->forNeighborsOf(v, [&](NetworKit::node u) {
        neighbors.push_back(u);
    });
    return neighbors;
}

std::vector<NetworKit::node> IndependentSet::getNeighborsPlus(NetworKit::node v) const {
    std::vector<NetworKit::node> neighborsPlus = getNeighbors(v);
    neighborsPlus.push_back(v);
    return neighborsPlus;
}

std::vector<NetworKit::node> IndependentSet::getNeighbors2(NetworKit::node v) const {    
    std::vector<NetworKit::node> neighbors = getNeighbors(v);
    std::vector<NetworKit::node> visitedNodes;    
    for (auto n : neighbors) {
        graph->forNeighborsOf(n, [&](NetworKit::node n2) {
            visitedNodes.push_back(n2);
        });
    }

    std::vector<bool> visitedOnGraph(graph->numberOfNodes());
    visitedOnGraph[v] = true;
    for (auto x : neighbors) {
        visitedOnGraph[x] = true;
    }

    std::vector<NetworKit::node> neighbors2;
    for (auto x : visitedNodes) {
        if (!visitedOnGraph[x]) {
            visitedOnGraph[x] = true;
            neighbors2.push_back(x);
        }
    }
    return neighbors2;
}

IndependentSet::EdgeSet IndependentSet::getConnectedEdges(std::vector<NetworKit::node>& nodes) const {
    EdgeSet connectedEdges(edgeComparator);
    for (auto x : nodes) {
        graph->forEdgesOf(x, [&](NetworKit::node u, NetworKit::node v) {
            connectedEdges.insert(NetworKit::Edge(u, v, true));
        });
    }
    return connectedEdges;
}

IndependentSet::EdgeSet IndependentSet::getInducedEdges(std::vector<NetworKit::node>& nodes) const {
    EdgeSet connectedEdges(edgeComparator); 
    std::set<NetworKit::node> nodeSet(nodes.begin(), nodes.end());

    for (auto u : nodes) {
        graph->forEdgesOf(u, [&](NetworKit::node u, NetworKit::node v) {
            if (nodeSet.contains(v)) {
                connectedEdges.insert(NetworKit::Edge(u, v, true));
            }
        });
    }
    return connectedEdges;
}

IndependentSet::EdgeSet IndependentSet::getAllEdges() const {
    EdgeSet edges(edgeComparator);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        edges.insert(NetworKit::Edge(u, v, true));
    });
    return edges;
}

std::vector<NetworKit::node> IndependentSet::getAllNodes() const {
    std::vector<NetworKit::node> nodes;
    graph->forNodes([&](NetworKit::node v) {
        nodes.push_back(v);
    });
    return nodes;
}

std::vector<NetworKit::node> IndependentSet::getMirrors(NetworKit::node v) const {
    std::vector<NetworKit::node> mirrors;
    std::vector<NetworKit::node> neighbors2 = getNeighbors2(v);
    std::vector<NetworKit::node> neighborsV = getNeighbors(v);

    for (auto w : neighbors2) {
        std::vector<NetworKit::node> neighborsW = getNeighbors(w);
        std::set<NetworKit::node> potentialClique(neighborsV.begin(), neighborsV.end());
        for (auto node : neighborsW) {
            potentialClique.erase(node);
        }

        bool clique = true;
        for (auto a : potentialClique) {
            for (auto b : potentialClique) {
                if (a != b && !graph->hasEdge(a, b)) {
                    clique = false;
                }
            }
        }
        if (clique) {
            mirrors.push_back(w);
        }
    }
    return mirrors;
}

void IndependentSet::dfs(NetworKit::node v, std::vector<bool>& visited) {
    if (visited[v]) {
        return;
    }
    visited[v] = true;
    std::vector<NetworKit::node> neighbors = getNeighbors(v);
    for (auto n : neighbors) {
       dfs(n, visited);
    }
}

NetworKit::node IndependentSet::getMinimumDegreeNode() const {
    assert(!graph->isEmpty());
    NetworKit::count minDegree = std::numeric_limits<NetworKit::count>::max();
    NetworKit::node index;
    graph->forNodes([&](NetworKit::node v) {
        if (minDegree > graph->degree(v)) {
            minDegree = graph->degree(v);
            index = v;
        }
    });
    return index;
}

NetworKit::node IndependentSet::getMaximumDegreeNode() const {
    assert(!graph->isEmpty());
    NetworKit::count maxDegree = 0;
    NetworKit::node index = -1;
    graph->forNodes([&](NetworKit::node v) {
        if (maxDegree <= graph->degree(v)) {
            maxDegree = graph->degree(v);
            index = v;
        }
    });
    return index;
}

void IndependentSet::removeElements(std::vector<NetworKit::node> nodes)  {
    for (auto v : nodes) {
        graph->removeNode(v);
    }
}

void IndependentSet::restoreElements(
        std::vector<NetworKit::node>& nodes, 
        EdgeSet& edges) {
    for (auto v : nodes) {
        graph->restoreNode(v);
    }
    for (auto e : edges) {
        graph->addEdge(e.u, e.v);
    }
}

std::vector<NetworKit::node> IndependentSet::runIndependentSetDegree2() const {
    NetworKit::Graph graphDeg2(*graph);
    std::vector<NetworKit::node> independentSet;

    auto solvePath = [&] (NetworKit::node u) { // u -> v -> w
        while (graphDeg2.hasNode(u) && graphDeg2.degree(u) <= 1) {
            if (graphDeg2.degree(u) == 0) {
                independentSet.push_back(u);
                graphDeg2.removeNode(u);
                break;
            }

            independentSet.push_back(u);
            NetworKit::node v = graphDeg2.getIthNeighbor(u, 0);
            graphDeg2.removeNode(u);

            if (graphDeg2.degree(v) == 1) {
                NetworKit::node w = graphDeg2.getIthNeighbor(v, 0);
                u = w;
            }
            graphDeg2.removeNode(v);
        }
    };

    graphDeg2.forNodes([&](NetworKit::node u) {
        solvePath(u);
    });

    graphDeg2.forNodes([&](NetworKit::node u) { // t <- u -> v -> w
        if (graphDeg2.degree(u) == 2) {
            independentSet.push_back(u);
            NetworKit::node t = graphDeg2.getIthNeighbor(u, 0);
            graphDeg2.removeNode(t);
            NetworKit::node v = graphDeg2.getIthNeighbor(u, 0);
            graphDeg2.removeNode(u);

            if (graphDeg2.degree(v) > 0) {
                NetworKit::node w = graphDeg2.getIthNeighbor(v, 0);
                graphDeg2.removeNode(v);
                solvePath(w);
            }
            else {
                graphDeg2.removeNode(v);
            }            
        }
    });
    return independentSet;
}

IndependentSet::IndependentSet(const NetworKit::Graph &graph) : 
        graph(std::make_optional(graph)) { }

bool IndependentSet::edgeComparator(const NetworKit::Edge& a, const NetworKit::Edge& b) {
    return a.u < b.u || (a.u == b.u && a.v < b.v);
}

const std::map<NetworKit::node, bool>& IndependentSet::getIndependentSet() const {
    assureFinished();
    return independentSet;
}

void BruteForceIndependentSet::run() {
    if (graph->isEmpty()) {
        return;
    }

    std::vector<bool> testSet;
    graph->forNodes([&](NetworKit::node v) {
        testSet.push_back(false);
        independentSet[v] = false;
    });

    int best = 0;
    unsigned long long max = 1;
    max <<= graph->numberOfNodes();

    for (unsigned long long binary = 1; binary != max; ++binary) {
        std::bitset<8 * sizeof(unsigned long long)> testSet(binary);
        size_t testSetSize = testSet.count();
        if (testSetSize <= best) {
            continue;
        }

        bool illegalEdge = false;
        graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
            if(testSet[u] && testSet[v]) {
                illegalEdge = true;
                return;
            }
        });

        if (!illegalEdge) {
            for (int i = 0; i < testSet.size(); ++i) {
                independentSet[i] = testSet[i];
                best = testSetSize;
            }
        }
    }

    hasRun = true;
}

} /* namespace Koala */

/*
TODO:
make better functions for standard and mirror branching
mis2 has duplicated code for branching with v and without v and mirrors
branching through connected components doesn't improve polynomial complexity because of how NetworKit::Graph is implemented
add NodeSet besides EdgeSet
use erase_if (vector)
add dodyxgen docs

expected output:
time ./build/benchmark/benchmark_independent_set 2 < ~/Downloads/graph9.g6
SIZE 1: 1
SIZE 2: 1896
SIZE 3: 101267
SIZE 4: 142276
SIZE 5: 27107
SIZE 6: 1995
SIZE 7: 117
SIZE 8: 8
SIZE 9: 1

real	0m2,069s
user	0m2,069s
sys	0m0,000s

*/
