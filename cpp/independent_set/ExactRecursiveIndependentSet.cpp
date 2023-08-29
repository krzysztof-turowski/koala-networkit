/*
 * ExactRecrusvieIndependentSet.cpp
 *
 *  Created on: 26.08.2023
 *      Author: Artur Salawa
 */

#include <independent_set/IndependentSet.hpp>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <limits>

namespace Koala {

void RecursiveIndependentSet::run() {
    std::vector<NetworKit::node> result = recursive();
    independentSet = std::set(result.begin(), result.end());
    hasRun = true;
}

std::vector<NetworKit::node> Mis1IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }
    NetworKit::node v = getMinimumDegreeNode();

    int selectedToSet;
    std::vector<NetworKit::node> largestSet;
    for (auto u : getNeighborsPlus(v)) {
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
            std::vector<NetworKit::node> neighbors2 = getNeighbors2(v);
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
        std::vector<NetworKit::Edge> allEdges(graph->edgeRange().begin(), graph->edgeRange().end());

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

    NetworKit::node v = getMaximumDegreeNode();
    if (graph->degree(v) >= 3) {
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

std::vector<NetworKit::node> Mis4IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }

    NetworKit::node v = getMaximumDegreeNode();
    if (graph->degree(v) >= 3) {
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

std::vector<NetworKit::node> Mis5IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }

    NetworKit::node v = getMaximumDegreeNode();
    if (graph->degree(v) >= 3) {
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

std::vector<NetworKit::node> MeasureAndConquerIndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }

    NetworKit::node lastIndex = std::numeric_limits<NetworKit::node>::min();
    graph->forNodes([&](NetworKit::node u) {
        lastIndex = u;
    });
    std::vector<bool> visited(lastIndex + 1, false); 
    dfs(lastIndex, visited);

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
        std::vector<NetworKit::Edge> allEdges(graph->edgeRange().begin(), graph->edgeRange().end());

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

    struct NodeNeighbors {
        NetworKit::node index;
        std::set<NetworKit::node> neighbors;
    };
    std::vector<NodeNeighbors> nodeNeighbors;
    graph->forNodes([&](NetworKit::node v) {
        std::vector<NetworKit::node> vec = getNeighborsPlus(v);
        nodeNeighbors.push_back({v, std::set<NetworKit::node>(vec.begin(), vec.end())});        
    });

    for (auto& nodeV: nodeNeighbors) {
        for (auto& nodeW: nodeNeighbors) {
            if (nodeV.index != nodeW.index) {
                if (std::includes(nodeV.neighbors.begin(), nodeV.neighbors.end(), 
                                  nodeW.neighbors.begin(), nodeW.neighbors.end())) { 
                    std::vector<NetworKit::node> justV = {nodeV.index};
                    EdgeSet connectedEdges = getConnectedEdges(justV);
                    removeElements(justV);
                    std::vector<NetworKit::node> setWithoutV = recursive();
                    restoreElements(justV, connectedEdges);
                    return setWithoutV;
                }
            }
        }
    } // vertices degree 1 don't exist at this point

    // A node v is foldable if N(v) contains no anti-triangles.
    // If none of the conditions above holds, and (3) there is a foldable node v of degree
    // d(v) <= 3, or a foldable node of degree d(v) = 4 with at most three anti-edges in
    // N(v), the algorithm selects one such node v of minimum degree and folds it
    for (auto v : graph->nodeRange()) {
        if (graph->degree(v) <= 4) {
            std::vector<NetworKit::node> vNeighbors = getNeighbors(v);
            EdgeSet inducedEdges = getInducedEdges(vNeighbors);
            bool shallFold = false;

            switch(graph->degree(v)) {
            case 0: {
                graph->removeNode(v);
                std::vector<NetworKit::node> setWithV = recursive();
                setWithV.push_back(v);
                graph->restoreNode(v);
                return setWithV;
                break;
            }
            case 1: {
                std::cout << "CAN'T BE HERE !!! (a)" << std::endl;
                break;
            }
            case 2: {
                //shallFold = true; break; // TODO also works
                if (inducedEdges.empty()) {
                    NetworKit::node u1 = vNeighbors[0];
                    NetworKit::node u2 = vNeighbors[1];
                    NetworKit::node u12 = u1;                    
                    std::vector<NetworKit::node> vNeighbors2 = getNeighbors2(v);

                    std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
                    EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
                    removeElements(neighborsPlus);
                    graph->restoreNode(u12);
                    for (auto u : vNeighbors2) {
                        graph->addEdge(u12, u);
                    }

                    std::vector<NetworKit::node> resultSet = recursive();
                    bool newNodeChoosen = false;
                    for (int i = 0; i < resultSet.size(); ++i) {
                        if (resultSet[i] == u12) {
                            resultSet.push_back(u2);
                            newNodeChoosen = true;
                            break;
                        }
                    }
                    if (!newNodeChoosen) {
                        resultSet.push_back(v);
                    }

                    graph->removeNode(u12);
                    restoreElements(neighborsPlus, connectedEdges);
                    return resultSet;
                }
                else {
                    std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
                    EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
                    removeElements(neighborsPlus);
                    std::vector<NetworKit::node> setWithV = recursive();
                    setWithV.push_back(v);
                    restoreElements(neighborsPlus, connectedEdges);
                    return setWithV;
                }
                break;
            }
            case 3: {
                if (!inducedEdges.empty()) {
                    if (inducedEdges.size() == 1) { // other cases removed by domination
                        shallFold = true;   
                    }
                    else {
                        std::cout << "CAN'T BE HERE !!! (b)" << std::endl;
                    }
                }
                break;
            }
            case 4: {
                if (6 - inducedEdges.size() <= 3) { // at most 3 anti-edges from K4
                    break;
                }

                if (inducedEdges.size() >= 4) { // no triangle possible with 4 edges
                    std::cout << "induced size >= 4 case reached" << std::endl; // TODO: tests never reach this
                    shallFold = true;
                }
                else if (inducedEdges.size() == 3) {
                    std::cout << "induced size == 3 case reached" << std::endl; // TODO: tests never reach this
                    // FOLD <=> no induced vertex has degree 3
                    bool hasDegree3 = false;
                    std::map<NetworKit::node, int> counter;                    
                    for (auto e : inducedEdges) {
                        ++counter[e.u];
                        ++counter[e.v];
                        if (counter[e.u] == 3 || counter[e.v] == 3) {
                            hasDegree3 = true;
                        }
                    }
                    shallFold = !hasDegree3;
                }
                break;
            }
            }

            if (shallFold) {
                struct FoldingNode {
                    NetworKit::node oldU1;
                    NetworKit::node oldU2;
                    NetworKit::node u12;                            
                };
                std::vector<FoldingNode> foldingNodes; // create one for each antiedge
                std::vector<std::set<NetworKit::node>> vNeighborsNeighbors(lastIndex + 1);
                std::set<NetworKit::node> vNeighborsPlusSet(vNeighbors.begin(), vNeighbors.end());
                vNeighborsPlusSet.insert(v);

                for (auto u : vNeighbors) {
                    std::vector<NetworKit::node> uNeighborsVect = getNeighbors(u);
                    std::set<NetworKit::node> uNeighborsSet(uNeighborsVect.begin(), uNeighborsVect.end());
                    std::set_difference(uNeighborsSet.begin(), uNeighborsSet.end(),
                                        vNeighborsPlusSet.begin(), vNeighborsPlusSet.end(), 
                                        std::inserter(vNeighborsNeighbors[u], vNeighborsNeighbors[u].end()));
                }              

                for (auto u : vNeighbors) {
                    for (auto v : vNeighbors) {
                        if (u < v) {
                            if (!graph->hasEdge(u, v)) {
                                foldingNodes.push_back({u, v});
                            }
                        }
                    }
                }
                for (int i = 0; i < foldingNodes.size(); ++i) {
                    foldingNodes[i].u12 = vNeighbors[i];
                }

                std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
                EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
                removeElements(neighborsPlus);

                for (auto f : foldingNodes) {
                    graph->restoreNode(f.u12);
                    for (auto node : vNeighborsNeighbors[f.oldU1]) {
                        graph->addEdge(f.u12, node);
                    }
                    for (auto node : vNeighborsNeighbors[f.oldU2]) {
                        if (!graph->hasEdge(f.u12, node)) {
                            graph->addEdge(f.u12, node);
                        }
                    }
                }
                for (auto f : foldingNodes) {
                    for (auto g : foldingNodes) {
                        if (f.u12 < g.u12) {
                            graph->addEdge(f.u12, g.u12);
                        }
                    }
                }

                std::vector<NetworKit::node> resultSet = recursive();
                bool newNodeChoosen = false;
                for (int i = 0; !newNodeChoosen && i < resultSet.size(); ++i) {
                    for (auto f : foldingNodes) {
                        if (resultSet[i] == f.u12) {
                            resultSet[i] = f.oldU1;
                            resultSet.push_back(f.oldU2);
                            newNodeChoosen = true;
                            break;
                        }
                    }                           
                }
                if (!newNodeChoosen) {
                    resultSet.push_back(v);
                }
                for (auto f : foldingNodes) {
                    graph->removeNode(f.u12);
                }
                restoreElements(neighborsPlus, connectedEdges);
                return resultSet;
            }
        }
    }

    NetworKit::node v = getMaximumDegreeNode();
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