
/*
 * SimpleIndependentSet.cpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#include <independentSet/SimpleIndependentSet.hpp>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <limits>
#include <set>

namespace Koala {


void SimpleIndependentSet::getdependentNodes(
        NetworKit::node v, 
        std::vector<NetworKit::node>& dependentNodes) const {
    dependentNodes.push_back(v);
    graph->forNeighborsOf(v, [&](NetworKit::node u) {
        dependentNodes.push_back(u);
    });
}

void SimpleIndependentSet::getdependentElements(
        NetworKit::node v, 
        std::vector<NetworKit::node>& dependentNodes, 
        EdgeSet& dependentEdges) const {
    getdependentNodes(v, dependentNodes);
    for (auto x : dependentNodes) {
        graph->forEdgesOf(x, [&](NetworKit::node u, NetworKit::node v) {
            dependentEdges.insert(NetworKit::Edge(u, v, true));
        });
    }
}

NetworKit::count SimpleIndependentSet::getGraphsMaximumDegree() const {
    NetworKit::count maximumDegree = 0;
    graph->forNodes([&](NetworKit::node v) {
        maximumDegree = std::max(maximumDegree, graph->degree(v));
    });
    return maximumDegree;
}

NetworKit::node SimpleIndependentSet::getMinimumDegreeNode() const {
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

NetworKit::node SimpleIndependentSet::getMaximumDegreeNode() const {
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

void SimpleIndependentSet::removeElements(std::vector<NetworKit::node> nodes)  {
    for (auto v : nodes) {
        graph->removeNode(v);
    }
}

void SimpleIndependentSet::restoreElements(
        std::vector<NetworKit::node>& nodes, 
        EdgeSet& edges) {
    for (auto v : nodes) {
        graph->restoreNode(v);
    }
    for (auto e : edges) {
        graph->addEdge(e.u, e.v);
    }
}

std::vector<NetworKit::node> SimpleIndependentSet::runIndependentSetDegree2() const {
    NetworKit::Graph graphDeg2(graph.value());
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
        if (graphDeg2.hasNode(u) && graphDeg2.degree(u) == 2) {
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

SimpleIndependentSet::SimpleIndependentSet(const NetworKit::Graph &graph) 
    : graph(std::make_optional(graph))
    , edgeComparator([](NetworKit::Edge a, NetworKit::Edge b) {
        return a.u < b.u || (a.u == b.u && a.v < b.v);
    }) { }

const std::map<NetworKit::node, bool>& SimpleIndependentSet::getIndependentSet() const {
    assureFinished();
    return independentSet;
}

void BruteForceIndependentSet::run() {
    if (graph->numberOfNodes() == 0) {
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
        unsigned long long tmp = binary;
        int testSetSize = 0;
        for (int i = 0; i < testSet.size(); ++i) {
            testSet[i] = (tmp % 2);
            testSetSize += (tmp % 2);
            tmp >>= 1;
        }

        if (testSetSize <= best) {
            continue;
        }

        bool bad = false;
        graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
            if(testSet[u] && testSet[v]) {
                bad = true;
                return;
            }
        });

        if (!bad) {
            for (int i = 0; i < testSet.size(); ++i) {
                independentSet[i] = testSet[i];
                best = testSetSize;
            }
        }
    }
    hasRun = true;
}

void Mis1IndependentSet::run() {
    std::vector<NetworKit::node> result = recursive();
    graph->forNodes([&](NetworKit::node v) {
        independentSet[v] = false;
    });
    for (auto node : result) {
        independentSet[node] = true;
    }   
    hasRun = true;
}
    
std::vector<NetworKit::node> Mis1IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }
    NetworKit::node v = getMinimumDegreeNode();
    std::vector<NetworKit::node> selectOneOf;
    getdependentNodes(v, selectOneOf);

    int selectedToSet;
    std::vector<NetworKit::node> largestSet;
    for (auto u : selectOneOf) {
        std::vector<NetworKit::node> dependentNodes; 
        EdgeSet dependentEdges(edgeComparator);
        getdependentElements(u, dependentNodes, dependentEdges);
        
        removeElements(dependentNodes);
        std::vector<NetworKit::node> bestWithU = recursive();
        if (bestWithU.size() >= largestSet.size()) {
            largestSet = bestWithU;
            selectedToSet = u;
        }
        restoreElements(dependentNodes, dependentEdges);
    }

    largestSet.push_back(selectedToSet);
    return largestSet;
}

void Mis2IndependentSet::run() {
    std::vector<NetworKit::node> result = recursive();
    graph->forNodes([&](NetworKit::node v) {
        independentSet[v] = false;
    });
    for (auto node : result) { 
        independentSet[node] = true;
    }   
    hasRun = true;
}

std::vector<NetworKit::node> Mis2IndependentSet::recursive() {
    return {};
}

void Mis3IndependentSet::run() {
    std::vector<NetworKit::node> result = recursive();
    graph->forNodes([&](NetworKit::node v) {
        independentSet[v] = false;
    });
    for (auto node : result) {
        independentSet[node] = true;
    }   
    hasRun = true;
}

std::vector<NetworKit::node> Mis3IndependentSet::recursive() {
    return {};
}

void Mis4IndependentSet::run() {
    std::vector<NetworKit::node> result = recursive();
    graph->forNodes([&](NetworKit::node v) {
        independentSet[v] = false;
    });
    for (auto node : result) {
        independentSet[node] = true;
    }   
    hasRun = true;
}

std::vector<NetworKit::node> Mis4IndependentSet::recursive() {
    return {};
}

void Mis5IndependentSet::run() {
    std::vector<NetworKit::node> result = recursive();
    graph->forNodes([&](NetworKit::node v) {
        independentSet[v] = false;
    });
    for (auto node : result) {
        independentSet[node] = true;
    }
    hasRun = true;
}

std::vector<NetworKit::node> Mis5IndependentSet::recursive() {
    if (getGraphsMaximumDegree() >= 3) {
        NetworKit::node v = getMaximumDegreeNode();

        std::vector<NetworKit::node> dependentNodes; 
        EdgeSet dependentEdges(edgeComparator);
        getdependentElements(v, dependentNodes, dependentEdges);
        removeElements(dependentNodes);
        std::vector<NetworKit::node> setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(dependentNodes, dependentEdges);
        
        std::vector<NetworKit::node> vNeighbors;
        graph->forEdgesOf(v, [&](NetworKit::node v, NetworKit::node u) {
            vNeighbors.push_back(u);
        });
        graph->removeNode(v);
        std::vector<NetworKit::node> setWithoutV = recursive();
        graph->restoreNode(v);
        for (auto u : vNeighbors) {
            graph->addEdge(v, u);
        }
        
        return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
    }
    else {
        return runIndependentSetDegree2();
    }
}

} /* namespace Koala */
