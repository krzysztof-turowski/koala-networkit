/*
 * IndependentSet.cpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#include <independent_set/IndependentSet.hpp>

#include <bitset>

namespace Koala {

IndependentSet::IndependentSet(const NetworKit::Graph &graph)
    : graph(std::make_optional(graph)) { }

bool IndependentSet::edgeComparator(const NetworKit::Edge& a, const NetworKit::Edge& b) {
    return a.u < b.u || (a.u == b.u && a.v < b.v);
}

const std::set<NetworKit::node>& IndependentSet::getIndependentSet() const {
    assureFinished();
    return independentSet;
}

void IndependentSet::check() const {
    assureFinished();
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        assert(!(independentSet.contains(u) && independentSet.contains(v)));
    });
}

std::vector<NetworKit::node> IndependentSet::getNeighbors(NetworKit::node v) const {
    return std::vector<NetworKit::node>(
        graph->neighborRange(v).begin(), graph->neighborRange(v).end());
}

std::vector<NetworKit::node> IndependentSet::getNeighborsPlus(NetworKit::node v) const {
    std::vector<NetworKit::node> neighborsPlus = getNeighbors(v);
    neighborsPlus.push_back(v);
    return neighborsPlus;
}

std::vector<NetworKit::node> IndependentSet::getNeighbors2(NetworKit::node v) const {
    std::vector<NetworKit::node> neighbors = getNeighbors(v), visitedNodes;
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

IndependentSet::EdgeSet IndependentSet::getConnectedEdges(
        std::vector<NetworKit::node>& nodes) const {
    EdgeSet connectedEdges(edgeComparator);
    for (auto u : nodes) {
        graph->forNeighborsOf(u, [&](NetworKit::node v) {
            connectedEdges.insert(NetworKit::Edge(u, v, true));
        });
    }
    return connectedEdges;
}

IndependentSet::EdgeSet IndependentSet::getInducedEdges(
        std::vector<NetworKit::node>& nodes) const {
    EdgeSet connectedEdges(edgeComparator);
    std::set<NetworKit::node> nodeSet(nodes.begin(), nodes.end());
    for (auto u : nodes) {
        graph->forNeighborsOf(u, [&](NetworKit::node v) {
            if (nodeSet.contains(v)) {
                connectedEdges.insert(NetworKit::Edge(u, v, true));
            }
        });
    }
    return connectedEdges;
}

std::vector<NetworKit::node> IndependentSet::getMirrors(NetworKit::node v) const {
    std::vector<NetworKit::node> mirrors, neighborsV = getNeighbors(v);
    for (auto w : getNeighbors2(v)) {
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

NetworKit::node IndependentSet::getMinimumDegreeNode() const {
    assert(!graph->isEmpty());
    return *std::min_element(graph->nodeRange().begin(), graph->nodeRange().end(),
        [&](NetworKit::node v, NetworKit::node u) {return graph->degree(v) < graph->degree(u);});
}

NetworKit::node IndependentSet::getMaximumDegreeNode() const {
    assert(!graph->isEmpty());
    return *std::max_element(graph->nodeRange().begin(), graph->nodeRange().end(),
        [&](NetworKit::node v, NetworKit::node u) {return graph->degree(v) < graph->degree(u);});
}

void IndependentSet::removeElements(std::vector<NetworKit::node> nodes) {
    for (auto v : nodes) {
        graph->removeNode(v);
    }
}

std::vector<NetworKit::node> IndependentSet::runIndependentSetDegree2() const {
    NetworKit::Graph graphDeg2(*graph);
    std::vector<NetworKit::node> independentSet;
    auto solvePath = [&](NetworKit::node u) {  // u -> v -> w
        while (graphDeg2.hasNode(u) && graphDeg2.degree(u) <= 1) {
            if (graphDeg2.degree(u) == 0) {
                independentSet.push_back(u);
                graphDeg2.removeNode(u);
                break;
            }
            independentSet.push_back(u);
            auto v = graphDeg2.getIthNeighbor(u, 0);
            graphDeg2.removeNode(u);
            if (graphDeg2.degree(v) == 1) {
                auto w = graphDeg2.getIthNeighbor(v, 0);
                u = w;
            }
            graphDeg2.removeNode(v);
        }
    };
    graphDeg2.forNodes([&](NetworKit::node u) {
        solvePath(u);
    });
    graphDeg2.forNodes([&](NetworKit::node u) {  // t <- u -> v -> w
        if (graphDeg2.degree(u) == 2) {
            independentSet.push_back(u);
            auto t = graphDeg2.getIthNeighbor(u, 0);
            graphDeg2.removeNode(t);
            auto v = graphDeg2.getIthNeighbor(u, 0);
            graphDeg2.removeNode(u);
            if (graphDeg2.degree(v) > 0) {
                auto w = graphDeg2.getIthNeighbor(v, 0);
                graphDeg2.removeNode(v);
                solvePath(w);
            } else {
                graphDeg2.removeNode(v);
            }
        }
    });
    return independentSet;
}

void BruteForceIndependentSet::run() {
    if (graph->isEmpty()) {
        return;
    }
    std::vector<bool> testSet;
    graph->forNodes([&](NetworKit::node v) {
        testSet.push_back(false);
    });
    int best = 0;
    uint64_t max = (1 << graph->numberOfNodes());
    for (unsigned long long binary = 1; binary != max; ++binary) {
        std::bitset<8 * sizeof(uint64_t)> testSet(binary);
        size_t testSetSize = testSet.count();
        if (testSetSize <= best) {
            continue;
        }
        bool illegalEdge = false;
        for (const auto& [u, v] : graph->edgeRange()) {
            if (testSet[u] && testSet[v]) {
                illegalEdge = true;
                break;
            }
        }
        if (!illegalEdge) {
            independentSet.clear();
            for (int i = 0; i < testSet.size(); ++i) {
                if (testSet[i]) {
                    independentSet.insert(i);
                }
            }
            best = testSetSize;
        }
    }
    hasRun = true;
}

}  /* namespace Koala */
