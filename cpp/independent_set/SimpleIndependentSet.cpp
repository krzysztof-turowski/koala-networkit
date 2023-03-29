
/*
 * SimpleIndependentSet.cpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#include <independent_set/SimpleIndependentSet.hpp>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <limits>
#include <set>

namespace Koala {

std::vector<NetworKit::node> SimpleIndependentSet::getNeighbors(NetworKit::node v) const {
    std::vector<NetworKit::node> neighborsPlus;
    graph->forNeighborsOf(v, [&](NetworKit::node u) {
        neighborsPlus.push_back(u);
    });
    return neighborsPlus;
}

std::vector<NetworKit::node> SimpleIndependentSet::getNeighborsPlus(NetworKit::node v) const {
    std::vector<NetworKit::node> neighborsPlus;
    neighborsPlus.push_back(v);
    graph->forNeighborsOf(v, [&](NetworKit::node u) {
        neighborsPlus.push_back(u);
    });
    return neighborsPlus;
}

std::set<NetworKit::node> SimpleIndependentSet::getNeighbors2(NetworKit::node v) const {
    std::set<NetworKit::node> neighbors2;
    std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
    for (auto x : neighborsPlus) {
        graph->forEdgesOf(x, [&](NetworKit::node u, NetworKit::node v) {
            neighbors2.insert(v);
        });
    }
    for (auto x : neighborsPlus) {
        neighbors2.erase(x);
    }
    return neighbors2;
}

std::set<NetworKit::node> SimpleIndependentSet::getNeighbors2Plus(NetworKit::node v) const {
    std::set<NetworKit::node> neighbors2Plus;
    std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);        
    for (auto x : neighborsPlus) {
        graph->forEdgesOf(x, [&](NetworKit::node u, NetworKit::node v) {
            neighbors2Plus.insert(v);
        });
    }
    return neighbors2Plus;
}

SimpleIndependentSet::EdgeSet SimpleIndependentSet::getConnectedEdges(std::vector<NetworKit::node>& nodes) {
    EdgeSet connectedEdges(edgeComparator);
    for (auto x : nodes) {
        graph->forEdgesOf(x, [&](NetworKit::node u, NetworKit::node v) {
            connectedEdges.insert(NetworKit::Edge(u, v, true));
        });
    }
    return connectedEdges;
}

std::vector<NetworKit::node> SimpleIndependentSet::getMirrors(NetworKit::node v) const {
    std::vector<NetworKit::node> mirrors;
    std::set<NetworKit::node> neighbors2 = getNeighbors2(v);
    std::vector<NetworKit::node> neighborsV = getNeighbors(v);

    for (auto w : neighbors2) {
        std::vector<NetworKit::node> neighborsW = getNeighbors(v);
        std::set<NetworKit::node> potentialClique(neighborsV.begin(), neighborsV.end());
        for (auto node : neighborsW) {
            potentialClique.erase(node);
        }

        bool clique = true;
        for (auto a : potentialClique) { // TODO: this might be super slow
            for (auto b : potentialClique) {
                if (!graph->hasEdge(a, b)) {
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



NetworKit::count SimpleIndependentSet::getGraphsMaximumDegree() const {
    assert(!graph->isEmpty());
    NetworKit::count maximumDegree = 0;
    graph->forNodes([&](NetworKit::node v) {
        maximumDegree = std::max(maximumDegree, graph->degree(v));
    });
    return maximumDegree;
}

NetworKit::node SimpleIndependentSet::getMinimumDegreeNode() const {
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

NetworKit::node SimpleIndependentSet::getMaximumDegreeNode() const {
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
    std::vector<NetworKit::node> selectOneOf = getNeighborsPlus(v);

    int selectedToSet;
    std::vector<NetworKit::node> largestSet;
    for (auto u : selectOneOf) {
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(u);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
 
        removeElements(neighborsPlus);
        std::vector<NetworKit::node> bestWithU = recursive();
        if (bestWithU.size() >= largestSet.size()) {
            largestSet = bestWithU;
            selectedToSet = u;
        }
        restoreElements(neighborsPlus, connectedEdges);
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
    // if |graph| = 0 then
    if (graph->isEmpty()) {
        return {};
    }

    NetworKit::node v = getMinimumDegreeNode();
    switch(graph->degree(v)) {
    case 0:
    case 1: {
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

        removeElements(neighborsPlus);
        std::vector<NetworKit::node> independentSet = recursive();
        independentSet.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);        
        return independentSet;
        break;
    }
    case 2: {
        NetworKit::node u1 = graph->getIthNeighbor(v, 0);
        NetworKit::node u2 = graph->getIthNeighbor(v, 1);
        if (graph->hasEdge(u1, u2)) {
            std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
            EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

            removeElements(neighborsPlus);
            std::vector<NetworKit::node> independentSet = recursive();
            independentSet.push_back(v);
            restoreElements(neighborsPlus, connectedEdges);
            return independentSet;
        }
        else {
            std::set<NetworKit::node> neighbors2 = getNeighbors2(v);
            if (neighbors2.size() == 1) { // TODO: book adds second branch but in my opinion it is obsolete ?
                NetworKit::node w = *neighbors2.begin();
                std::vector<NetworKit::node> nodesToBeRemoved{v, u1, u2, w};
                EdgeSet connectedEdges = getConnectedEdges(nodesToBeRemoved);
                removeElements(nodesToBeRemoved);
                std::vector<NetworKit::node> independentSet = recursive();
                independentSet.push_back(u1);
                independentSet.push_back(u2);
                restoreElements(nodesToBeRemoved, connectedEdges);
                return independentSet;
            }
            else { // TODO: book doesn't say about +1 in the first case ?
                std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
                EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
                removeElements(neighborsPlus);
                std::vector<NetworKit::node> setWithTwoDegreeNode = recursive();
                setWithTwoDegreeNode.push_back(v);
                restoreElements(neighborsPlus, connectedEdges);

                std::vector<NetworKit::node> twoDegreeNodeWithMirrors = getMirrors(v);
                twoDegreeNodeWithMirrors.push_back(v);
                connectedEdges = getConnectedEdges(twoDegreeNodeWithMirrors);
                removeElements(twoDegreeNodeWithMirrors);
                std::vector<NetworKit::node> setWithoutTwoDegreeNode = recursive();
                setWithoutTwoDegreeNode.push_back(v);
                restoreElements(twoDegreeNodeWithMirrors, connectedEdges);

                return setWithTwoDegreeNode.size() > setWithoutTwoDegreeNode.size() ? 
                        setWithTwoDegreeNode : setWithoutTwoDegreeNode;
            }
        }
        break;
    }
    case 3:{

        break;
    }
    default: {

        break;
    }
    }

    // if exists v with d(v) >= 6 then

    //if graph is disconnected then

    //if graph is 4 or 5-regular then

    //if every v has d(v) \in {4, 5}

    // this is just a filler for all the not implemented cases
    if (getGraphsMaximumDegree() >= 3) {
        NetworKit::node v = getMaximumDegreeNode();
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

        removeElements(neighborsPlus);
        std::vector<NetworKit::node> setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);
        
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
    if (graph->isEmpty()) {
        return {};
    }

    std::optional<NetworKit::node> zeroDegreeNode;
    graph->forNodes([&](NetworKit::node v) {
        if (graph->degree(v) == 0) {
            zeroDegreeNode = v;
        }
    });
    if (zeroDegreeNode.has_value()) {
        graph->removeNode(*zeroDegreeNode);
        std::vector<NetworKit::node> independentSet = recursive();
        independentSet.push_back(*zeroDegreeNode);
        graph->restoreNode(*zeroDegreeNode);
        return independentSet;
    }

    std::optional<NetworKit::node> oneDegreeNode;
    graph->forNodes([&](NetworKit::node v) {
        if (graph->degree(v) == 1) {
            oneDegreeNode = v;
        }
    });
    if (oneDegreeNode.has_value()) {
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(*oneDegreeNode);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

        removeElements(neighborsPlus);
        std::vector<NetworKit::node> independentSet = recursive();
        independentSet.push_back(*oneDegreeNode);
        restoreElements(neighborsPlus, connectedEdges);        
        return independentSet;
    }

    if (getGraphsMaximumDegree() >= 3) {
        NetworKit::node v = getMaximumDegreeNode();
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        
        removeElements(neighborsPlus);
        std::vector<NetworKit::node> setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);
        
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
    if (graph->isEmpty()) {
        return {};
    }

    if (getGraphsMaximumDegree() >= 3) {
        NetworKit::node v;
        graph->forNodes([&](NetworKit::node u) {
            if (graph->degree(u) >= 3) {
                v = u;
            }
        });

        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

        removeElements(neighborsPlus);
        std::vector<NetworKit::node> setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);
        
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
    if (graph->isEmpty()) {
        return {};
    }

    if (getGraphsMaximumDegree() >= 3) {
        NetworKit::node v = getMaximumDegreeNode();

        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

        removeElements(neighborsPlus);
        std::vector<NetworKit::node> setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);
        
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

/*
TODO:
https://stackoverflow.com/questions/16491675/how-to-send-custom-message-in-google-c-testing-framework
maybe move run() function to the SimpleIndependentSet cuz its always the same and call overloaded recursive() functions
maybe make function for remove, recursion, restore 
maybe optimize neighbors2 functions, also mis2 is slow


add dodyxgen docs
*/