
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

void printGraph(const NetworKit::Graph* graph) { // TODO: its just debug
    graph->forNodes([&](NetworKit::node v) {
        std::cout << v << " ";
    });
    std::cout << std::endl;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        std::cout << "(" << u << "," << v << ") ";
    });
    std::cout << std::endl;
}

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

SimpleIndependentSet::EdgeSet SimpleIndependentSet::getConnectedEdges(std::vector<NetworKit::node>& nodes) const {
    EdgeSet connectedEdges(edgeComparator);
    for (auto x : nodes) {
        graph->forEdgesOf(x, [&](NetworKit::node u, NetworKit::node v) {
            connectedEdges.insert(NetworKit::Edge(u, v, true));
        });
    }
    return connectedEdges;
}

SimpleIndependentSet::EdgeSet SimpleIndependentSet::getAllEdges() const {
    EdgeSet edges(edgeComparator);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        edges.insert(NetworKit::Edge(u, v, true));
    });
    return edges;
}

std::vector<NetworKit::node> SimpleIndependentSet::getMirrors(NetworKit::node v) const {
    std::vector<NetworKit::node> mirrors;
    std::set<NetworKit::node> neighbors2 = getNeighbors2(v);
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

void SimpleIndependentSet::dfs(NetworKit::node v, std::vector<bool>& visited) {
    if (visited[v]) {
        return;
    }
    visited[v] = true;
    std::vector<NetworKit::node> neighbors = getNeighbors(v);
    for (auto n : neighbors) {
       dfs(n, visited);
    }
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
                std::vector<NetworKit::node> setWithV = recursive();
                setWithV.push_back(v);
                restoreElements(neighborsPlus, connectedEdges);

                std::vector<NetworKit::node> vWithMirrors = getMirrors(v);
                vWithMirrors.push_back(v);
                connectedEdges = getConnectedEdges(vWithMirrors);
                removeElements(vWithMirrors);
                std::vector<NetworKit::node> setWithoutV = recursive();
                restoreElements(vWithMirrors, connectedEdges);

                return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
            }
        }
        break;
    }
    case 3:{
        NetworKit::node u1 = graph->getIthNeighbor(v, 0);
        NetworKit::node u2 = graph->getIthNeighbor(v, 1);
        NetworKit::node u3 = graph->getIthNeighbor(v, 2);
        if (!graph->hasEdge(u1, u2) && !graph->hasEdge(u2, u3) && !graph->hasEdge(u1, u3)) {
            std::vector<NetworKit::node> vMirrors = getMirrors(v);
            if (!vMirrors.empty()) {
                std::vector<NetworKit::node> neighborsPlus = {v, u1, u2, u3};
                EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
                removeElements(neighborsPlus);
                std::vector<NetworKit::node> setWithV = recursive();
                setWithV.push_back(v);
                restoreElements(neighborsPlus, connectedEdges);

                std::vector<NetworKit::node> vWithMirrors = getMirrors(v);
                vWithMirrors.push_back(v);
                connectedEdges = getConnectedEdges(vWithMirrors);
                removeElements(vWithMirrors);
                std::vector<NetworKit::node> setWithoutV = recursive();
                restoreElements(vWithMirrors, connectedEdges);

                return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
            }
            else {
                std::vector<NetworKit::node> vNeighborsPlus = {v, u1, u2, u3};
                graph->removeNode(v);
                std::vector<NetworKit::node> u1NeighborsPlus = getNeighborsPlus(u1); // without v
                std::vector<NetworKit::node> u2NeighborsPlus = getNeighborsPlus(u2); // without v
                std::vector<NetworKit::node> u3NeighborsPlus = getNeighborsPlus(u3); // without v
                graph->restoreNode(v);
                graph->addEdge(v, u1);
                graph->addEdge(v, u2);
                graph->addEdge(v, u3);

                std::vector<NetworKit::node> delVerticesCase1 = vNeighborsPlus;
                std::vector<NetworKit::node> delVerticesCase2 = {v};
                delVerticesCase2.insert(delVerticesCase2.end(), u1NeighborsPlus.begin(), u1NeighborsPlus.end());
                delVerticesCase2.insert(delVerticesCase2.end(), u2NeighborsPlus.begin(), u2NeighborsPlus.end());
                std::vector<NetworKit::node> delVerticesCase3 = {v, u2};
                delVerticesCase3.insert(delVerticesCase3.end(), u1NeighborsPlus.begin(), u1NeighborsPlus.end());
                delVerticesCase3.insert(delVerticesCase3.end(), u3NeighborsPlus.begin(), u3NeighborsPlus.end());
                std::vector<NetworKit::node> delVerticesCase4 = {v, u1};
                delVerticesCase4.insert(delVerticesCase4.end(), u2NeighborsPlus.begin(), u2NeighborsPlus.end());
                delVerticesCase4.insert(delVerticesCase4.end(), u3NeighborsPlus.begin(), u3NeighborsPlus.end());

                EdgeSet connectedEdgesCase1 = getConnectedEdges(delVerticesCase1);
                EdgeSet connectedEdgesCase2 = getConnectedEdges(delVerticesCase2);
                EdgeSet connectedEdgesCase3 = getConnectedEdges(delVerticesCase3);
                EdgeSet connectedEdgesCase4 = getConnectedEdges(delVerticesCase4);

                std::vector<NetworKit::node> setCase1, setCase2, setCase3, setCase4;
                removeElements(delVerticesCase1);
                setCase1 = recursive();
                setCase1.push_back(v);
                restoreElements(delVerticesCase1, connectedEdgesCase1);

                removeElements(delVerticesCase2);
                setCase2 = recursive();
                setCase2.push_back(u1);
                setCase2.push_back(u2);
                restoreElements(delVerticesCase2, connectedEdgesCase2);
                
                removeElements(delVerticesCase3);                
                setCase3 = recursive();
                setCase3.push_back(u1);
                setCase3.push_back(u3);
                restoreElements(delVerticesCase3, connectedEdgesCase3);
                             
                removeElements(delVerticesCase4);
                setCase4 = recursive();
                setCase4.push_back(u2);
                setCase4.push_back(u3);
                restoreElements(delVerticesCase4, connectedEdgesCase4);
                
                if (setCase1.size() >= setCase2.size() && setCase1.size() >= setCase3.size() && setCase1.size() >= setCase4.size()) return setCase1;
                if (setCase2.size() >= setCase1.size() && setCase2.size() >= setCase3.size() && setCase2.size() >= setCase4.size()) return setCase2;
                if (setCase3.size() >= setCase1.size() && setCase3.size() >= setCase2.size() && setCase3.size() >= setCase4.size()) return setCase3;
                return setCase4;
            }
        }
        else if (graph->hasEdge(u1, u2) && graph->hasEdge(u2, u3) && graph->hasEdge(u1, u3)) {
            std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
            EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
            removeElements(neighborsPlus);
            std::vector<NetworKit::node> setWithV = recursive();
            setWithV.push_back(v);
            restoreElements(neighborsPlus, connectedEdges);
            return setWithV;
        }
        else { // G[{u1,u2,u3}] has 1 or 2 edges
            std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
            EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
            removeElements(neighborsPlus);
            std::vector<NetworKit::node> setWithV = recursive();
            setWithV.push_back(v);
            restoreElements(neighborsPlus, connectedEdges);

            std::vector<NetworKit::node> vWithMirrors = getMirrors(v);
            vWithMirrors.push_back(v);
            connectedEdges = getConnectedEdges(vWithMirrors);
            removeElements(vWithMirrors);
            std::vector<NetworKit::node> setWithoutV = recursive();
            restoreElements(vWithMirrors, connectedEdges);
            return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
        }

        break;
    }
    } // no verticle with deg 0,1,2,3 exists after leaving this switch

    v = getMaximumDegreeNode();
    if (graph->degree(v) >= 6) {
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        removeElements(neighborsPlus);
        std::vector<NetworKit::node> setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);
        
        std::vector<NetworKit::node> justV = {v};
        connectedEdges = getConnectedEdges(justV);
        removeElements(justV);
        std::vector<NetworKit::node> setWithoutV = recursive();
        restoreElements(justV, connectedEdges);
       
        return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
    }

    NetworKit::node lastIndex = std::numeric_limits<NetworKit::node>::min();
    graph->forNodes([&](NetworKit::node u) {
        lastIndex = u;
    });
    std::vector<bool> visited(lastIndex + 1, false); 
    dfs(v, visited);

    std::vector<NetworKit::node> component;
    std::vector<NetworKit::node> theRest;
    std::vector<NetworKit::node> allVertices;
    graph->forNodes([&](NetworKit::node u) {
        if (visited[u]) {
            component.push_back(u);
        }
        else {
            theRest.push_back(u);
        }
        allVertices.push_back(u);
    });
    
    if (component.size() != allVertices.size()) {
        EdgeSet allEdges = getAllEdges();

        removeElements(theRest);
        std::vector<NetworKit::node> setComponent = recursive();
        removeElements(component);
        restoreElements(allVertices, allEdges);

        removeElements(component);
        std::vector<NetworKit::node> setTheRest = recursive();
        removeElements(theRest);
        restoreElements(allVertices, allEdges);

        setComponent.insert(setComponent.end(), setTheRest.begin(), setTheRest.end());
        return setComponent;
    }
  
    bool hasDegree4 = false;
    bool hasDegree5 = false;
    graph->forNodes([&](NetworKit::node u) { // all must be 4 or 5
        if (graph->degree(u) == 4) {
            hasDegree4 = true;
        }
        else {
            hasDegree5 = true;
        }
    });

    if (!(hasDegree4 && hasDegree5)) {
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        removeElements(neighborsPlus);
        std::vector<NetworKit::node> setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);

        std::vector<NetworKit::node> vWithMirrors = getMirrors(v);
        vWithMirrors.push_back(v);
        connectedEdges = getConnectedEdges(vWithMirrors);
        removeElements(vWithMirrors);
        std::vector<NetworKit::node> setWithoutV = recursive();
        restoreElements(vWithMirrors, connectedEdges);

        return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
    }

    // TODO: book counts v twice in second case
    NetworKit::node w;
    graph->forEdges([&](NetworKit::node a, NetworKit::node b) {
        if (graph->degree(a) != graph->degree(b)) {
            if (graph->degree(a) == 5) {
                v = a;
                w = b;
            }
            else {
                v = b;
                w = a;
            }
        }        
    });

    std::vector<NetworKit::node> vNeighborsPlus = getNeighborsPlus(v); // for case 1

    std::vector<NetworKit::node> vMirrors = getMirrors(v); // for case 2
    std::vector<NetworKit::node> wNeighborsPlus = getNeighborsPlus(w);
    std::set<NetworKit::node> case2TemporarySet(vMirrors.begin(), vMirrors.end());
    case2TemporarySet.insert(wNeighborsPlus.begin(), wNeighborsPlus.end());
    std::vector<NetworKit::node> vMirrorsWNeighborsPlus(case2TemporarySet.begin(), case2TemporarySet.end());

    std::vector<NetworKit::node> vMirrorsPlusW = vMirrors; // for case 3
    vMirrorsPlusW.push_back(v);
    vMirrorsPlusW.push_back(w);

    EdgeSet connectedEdgesCase1 = getConnectedEdges(vNeighborsPlus);
    EdgeSet connectedEdgesCase2 = getConnectedEdges(vMirrorsWNeighborsPlus);
    EdgeSet connectedEdgesCase3 = getConnectedEdges(vMirrorsPlusW);

    std::vector<NetworKit::node> setCase1, setCase2, setCase3;
    removeElements(vNeighborsPlus);
    setCase1 = recursive();
    setCase1.push_back(v);
    restoreElements(vNeighborsPlus, connectedEdgesCase1);

    removeElements(vMirrorsWNeighborsPlus);
    setCase2 = recursive();
    setCase2.push_back(w);
    restoreElements(vMirrorsWNeighborsPlus, connectedEdgesCase2);
    
    removeElements(vMirrorsPlusW);                
    setCase3 = recursive();
    restoreElements(vMirrorsPlusW, connectedEdgesCase3);    

    if (setCase1.size() >= setCase2.size() && setCase1.size() >= setCase3.size()) return setCase1;
    if (setCase2.size() >= setCase1.size() && setCase2.size() >= setCase3.size()) return setCase2;
    return setCase3;
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
        
        std::vector<NetworKit::node> justV = {v};
        connectedEdges = getConnectedEdges(justV);
        removeElements(justV);
        std::vector<NetworKit::node> setWithoutV = recursive();
        restoreElements(justV, connectedEdges);
        
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
        
        std::vector<NetworKit::node> justV = {v};
        connectedEdges = getConnectedEdges(justV);
        removeElements(justV);
        std::vector<NetworKit::node> setWithoutV = recursive();
        restoreElements(justV, connectedEdges);
        
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
        
        std::vector<NetworKit::node> justV = {v};
        connectedEdges = getConnectedEdges(justV);
        removeElements(justV);
        std::vector<NetworKit::node> setWithoutV = recursive();
        restoreElements(justV, connectedEdges);
        
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
maybe optimize neighbors2 functions
rename AdjacencyMatrix in tests
mis2 has duplicated code for branching with v and without v and mirrors
branching through connected components doesn't improve polynomial complexity because of how NetworKit::Graph is implemented
max/small degree graph function exists in the library => use it

add dodyxgen docs
*/