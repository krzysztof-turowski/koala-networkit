/*
 * ExactRecrusvieIndependentSet.cpp
 *
 *  Created on: 26.08.2023
 *      Author: Artur Salawa
 */

#include <independent_set/ExactIndependentSet.hpp>

#include <networkit/graph/GraphTools.hpp>
#include <traversal/DFS.hpp>

#include <limits>

namespace Koala {

void RecursiveIndependentSet::run() {
    auto result = recursive();
    independentSet = std::set(result.begin(), result.end());
    hasRun = true;
}

std::vector<NetworKit::node> Mis1IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }
    int selectedToSet;
    std::vector<NetworKit::node> largestSet;
    for (auto u : getNeighborsPlus(getMinimumDegreeNode())) {
        auto neighborsPlus = getNeighborsPlus(u);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        removeElements(neighborsPlus);
        auto bestWithU = recursive();
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

    auto v = getMinimumDegreeNode();
    switch (graph->degree(v)) {
    case 0:
    case 1: {
        auto neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

        removeElements(neighborsPlus);
        auto independentSet = recursive();
        independentSet.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);
        return independentSet;
    }
    case 2: {
        auto u1 = graph->getIthNeighbor(v, 0), u2 = graph->getIthNeighbor(v, 1);
        if (graph->hasEdge(u1, u2)) {
            auto neighborsPlus = getNeighborsPlus(v);
            EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

            removeElements(neighborsPlus);
            auto independentSet = recursive();
            independentSet.push_back(v);
            restoreElements(neighborsPlus, connectedEdges);
            return independentSet;
        } else {
            auto neighbors2 = getNeighbors2(v);
            if (neighbors2.size() == 1) {  // book adds second branch but IMO it is obsolete
                NetworKit::node w = *neighbors2.begin();
                std::vector<NetworKit::node> nodesToBeRemoved{v, u1, u2, w};
                EdgeSet connectedEdges = getConnectedEdges(nodesToBeRemoved);
                removeElements(nodesToBeRemoved);
                auto independentSet = recursive();
                independentSet.push_back(u1), independentSet.push_back(u2);
                restoreElements(nodesToBeRemoved, connectedEdges);
                return independentSet;
            } else {
                auto neighborsPlus = getNeighborsPlus(v);
                EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
                removeElements(neighborsPlus);
                auto setWithV = recursive();
                setWithV.push_back(v);
                restoreElements(neighborsPlus, connectedEdges);
                auto vWithMirrors = getMirrors(v);
                vWithMirrors.push_back(v);
                connectedEdges = getConnectedEdges(vWithMirrors);
                removeElements(vWithMirrors);
                auto setWithoutV = recursive();
                restoreElements(vWithMirrors, connectedEdges);

                return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
            }
        }
        break;
    }
    case 3: {
        auto u1 = graph->getIthNeighbor(v, 0);
        auto u2 = graph->getIthNeighbor(v, 1);
        auto u3 = graph->getIthNeighbor(v, 2);
        if (!graph->hasEdge(u1, u2) && !graph->hasEdge(u2, u3) && !graph->hasEdge(u1, u3)) {
            std::vector<NetworKit::node> vMirrors = getMirrors(v);
            if (!vMirrors.empty()) {
                std::vector<NetworKit::node> neighborsPlus{v, u1, u2, u3};
                EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
                removeElements(neighborsPlus);
                auto setWithV = recursive();
                setWithV.push_back(v);
                restoreElements(neighborsPlus, connectedEdges);

                auto vWithMirrors = getMirrors(v);
                vWithMirrors.push_back(v);
                connectedEdges = getConnectedEdges(vWithMirrors);
                removeElements(vWithMirrors);
                auto setWithoutV = recursive();
                restoreElements(vWithMirrors, connectedEdges);

                return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
            } else {
                std::vector<NetworKit::node> vNeighborsPlus = {v, u1, u2, u3};
                graph->removeNode(v);
                auto u1NeighborsPlus = getNeighborsPlus(u1);  // without v
                auto u2NeighborsPlus = getNeighborsPlus(u2);  // without v
                auto u3NeighborsPlus = getNeighborsPlus(u3);  // without v
                graph->restoreNode(v);
                graph->addEdge(v, u1), graph->addEdge(v, u2), graph->addEdge(v, u3);

                auto delVerticesCase1 = vNeighborsPlus;
                auto delVerticesCase2 = std::vector<NetworKit::node>{v};
                delVerticesCase2.insert(
                    delVerticesCase2.end(), u1NeighborsPlus.begin(), u1NeighborsPlus.end());
                delVerticesCase2.insert(
                    delVerticesCase2.end(), u2NeighborsPlus.begin(), u2NeighborsPlus.end());
                std::vector<NetworKit::node> delVerticesCase3 = {v, u2};
                delVerticesCase3.insert(
                    delVerticesCase3.end(), u1NeighborsPlus.begin(), u1NeighborsPlus.end());
                delVerticesCase3.insert(
                    delVerticesCase3.end(), u3NeighborsPlus.begin(), u3NeighborsPlus.end());
                std::vector<NetworKit::node> delVerticesCase4 = {v, u1};
                delVerticesCase4.insert(
                    delVerticesCase4.end(), u2NeighborsPlus.begin(), u2NeighborsPlus.end());
                delVerticesCase4.insert(
                    delVerticesCase4.end(), u3NeighborsPlus.begin(), u3NeighborsPlus.end());

                EdgeSet connectedEdgesCase1 = getConnectedEdges(delVerticesCase1);
                EdgeSet connectedEdgesCase2 = getConnectedEdges(delVerticesCase2);
                EdgeSet connectedEdgesCase3 = getConnectedEdges(delVerticesCase3);
                EdgeSet connectedEdgesCase4 = getConnectedEdges(delVerticesCase4);

                removeElements(delVerticesCase1);
                auto setCase1 = recursive();
                setCase1.push_back(v);
                restoreElements(delVerticesCase1, connectedEdgesCase1);

                removeElements(delVerticesCase2);
                auto setCase2 = recursive();
                setCase2.push_back(u1), setCase2.push_back(u2);
                restoreElements(delVerticesCase2, connectedEdgesCase2);

                removeElements(delVerticesCase3);
                auto setCase3 = recursive();
                setCase3.push_back(u1), setCase3.push_back(u3);
                restoreElements(delVerticesCase3, connectedEdgesCase3);

                removeElements(delVerticesCase4);
                auto setCase4 = recursive();
                setCase4.push_back(u2), setCase4.push_back(u3);
                restoreElements(delVerticesCase4, connectedEdgesCase4);

                std::vector<std::vector<NetworKit::node>> solutions{
                    setCase1, setCase2, setCase3, setCase4};
                return *std::max_element(
                    solutions.begin(), solutions.end(),
                    [](std::vector<NetworKit::node> a, std::vector<NetworKit::node> b)
                    {return a.size() < b.size();});
            }
        } else if (graph->hasEdge(u1, u2) && graph->hasEdge(u2, u3) && graph->hasEdge(u1, u3)) {
            auto neighborsPlus = getNeighborsPlus(v);
            EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
            removeElements(neighborsPlus);
            std::vector<NetworKit::node> setWithV = recursive();
            setWithV.push_back(v);
            restoreElements(neighborsPlus, connectedEdges);
            return setWithV;
        } else {  // G[{u1,u2,u3}] has 1 or 2 edges
            auto neighborsPlus = getNeighborsPlus(v);
            EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
            removeElements(neighborsPlus);
            auto setWithV = recursive();
            setWithV.push_back(v);
            restoreElements(neighborsPlus, connectedEdges);
            auto vWithMirrors = getMirrors(v);
            vWithMirrors.push_back(v);
            connectedEdges = getConnectedEdges(vWithMirrors);
            removeElements(vWithMirrors);
            auto setWithoutV = recursive();
            restoreElements(vWithMirrors, connectedEdges);
            return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
        }
        break;
    }
    }  // no verticle with deg 0,1,2,3 exists after leaving this switch

    v = getMaximumDegreeNode();
    if (graph->degree(v) >= 6) {
        auto neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        removeElements(neighborsPlus);
        auto setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);
        auto justV = std::vector<NetworKit::node>{v};
        connectedEdges = getConnectedEdges(justV);
        removeElements(justV);
        auto setWithoutV = recursive();
        restoreElements(justV, connectedEdges);

        return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
    }

    std::vector<bool> visited(graph->upperNodeIdBound(), false);
    Koala::Traversal::DFSFrom(
        *graph, *graph->nodeRange().begin(),
        [&](auto v) { visited[v] = true; }, [&](auto v) { return true; });
    std::vector<NetworKit::node> component, theRest, allVertices;
    graph->forNodes([&](NetworKit::node u) {
        if (visited[u]) {
            component.push_back(u);
        } else {
            theRest.push_back(u);
        }
        allVertices.push_back(u);
    });

    if (component.size() != allVertices.size()) {
        std::vector<NetworKit::Edge> allEdges(graph->edgeRange().begin(), graph->edgeRange().end());

        removeElements(theRest);
        auto setComponent = recursive();
        removeElements(component);
        restoreElements(allVertices, allEdges);

        removeElements(component);
        auto setTheRest = recursive();
        removeElements(theRest);
        restoreElements(allVertices, allEdges);

        setComponent.insert(setComponent.end(), setTheRest.begin(), setTheRest.end());
        return setComponent;
    }

    bool hasDegree4 = false, hasDegree5 = false;
    graph->forNodes([&](NetworKit::node u) {  // all must be 4 or 5
        if (graph->degree(u) == 4) {
            hasDegree4 = true;
        } else {
            hasDegree5 = true;
        }
    });

    if (!(hasDegree4 && hasDegree5)) {
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        removeElements(neighborsPlus);
        auto setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);

        auto vWithMirrors = getMirrors(v);
        vWithMirrors.push_back(v);
        connectedEdges = getConnectedEdges(vWithMirrors);
        removeElements(vWithMirrors);
        auto setWithoutV = recursive();
        restoreElements(vWithMirrors, connectedEdges);

        return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
    }
    NetworKit::node w;
    for (const auto &[a, b] : graph->edgeRange()) {
        if (graph->degree(a) != graph->degree(b)) {
            if (graph->degree(a) == 5) {
                v = a, w = b;
            } else {
                v = b, w = a;
            }
            break;
        }
    }

    auto vNeighborsPlus = getNeighborsPlus(v);  // for case 1
    auto vMirrors = getMirrors(v);  // for case 2
    auto wNeighborsPlus = getNeighborsPlus(w);
    std::set<NetworKit::node> case2TemporarySet(vMirrors.begin(), vMirrors.end());
    case2TemporarySet.insert(wNeighborsPlus.begin(), wNeighborsPlus.end());
    std::vector<NetworKit::node> vMirrorsWNeighborsPlus(
        case2TemporarySet.begin(), case2TemporarySet.end());

    auto vMirrorsPlusW = vMirrors;  // for case 3
    vMirrorsPlusW.push_back(v), vMirrorsPlusW.push_back(w);

    EdgeSet connectedEdgesCase1 = getConnectedEdges(vNeighborsPlus);
    EdgeSet connectedEdgesCase2 = getConnectedEdges(vMirrorsWNeighborsPlus);
    EdgeSet connectedEdgesCase3 = getConnectedEdges(vMirrorsPlusW);

    removeElements(vNeighborsPlus);
    auto setCase1 = recursive();
    setCase1.push_back(v);
    restoreElements(vNeighborsPlus, connectedEdgesCase1);

    removeElements(vMirrorsWNeighborsPlus);
    auto setCase2 = recursive();
    setCase2.push_back(w);
    restoreElements(vMirrorsWNeighborsPlus, connectedEdgesCase2);

    removeElements(vMirrorsPlusW);
    auto setCase3 = recursive();
    restoreElements(vMirrorsPlusW, connectedEdgesCase3);

    if (setCase1.size() >= setCase2.size() && setCase1.size() >= setCase3.size()) {
        return setCase1;
    }
    if (setCase2.size() >= setCase1.size() && setCase2.size() >= setCase3.size()) {
        return setCase2;
    }
    return setCase3;
}

std::vector<NetworKit::node> Mis3IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }
    auto u = getMinimumDegreeNode();
    switch (graph->degree(u)) {
    case 0: {
        graph->removeNode(u);
        std::vector<NetworKit::node> independentSet = recursive();
        independentSet.push_back(u);
        graph->restoreNode(u);
        return independentSet;
    }
    case 1: {
        std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(u);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);

        removeElements(neighborsPlus);
        auto independentSet = recursive();
        independentSet.push_back(u);
        restoreElements(neighborsPlus, connectedEdges);
        return independentSet;
    }
    default:
        auto v = getMaximumDegreeNode();
        if (graph->degree(v) >= 3) {
            auto v = getMaximumDegreeNode();
            std::vector<NetworKit::node> neighborsPlus = getNeighborsPlus(v);
            EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
            removeElements(neighborsPlus);
            auto setWithV = recursive();
            setWithV.push_back(v);
            restoreElements(neighborsPlus, connectedEdges);

            auto justV = std::vector<NetworKit::node>{v};
            connectedEdges = getConnectedEdges(justV);
            removeElements(justV);
            auto setWithoutV = recursive();
            restoreElements(justV, connectedEdges);

            return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
        } else {
            return runIndependentSetDegree2();
        }
    }
}

std::vector<NetworKit::node> Mis4IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }

    auto v = getMaximumDegreeNode();
    if (graph->degree(v) >= 3) {
        graph->forNodes([&](NetworKit::node u) {
            if (graph->degree(u) >= 3) {
                v = u;
            }
        });

        auto neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        removeElements(neighborsPlus);
        auto setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);

        auto justV = std::vector<NetworKit::node>{v};
        connectedEdges = getConnectedEdges(justV);
        removeElements(justV);
        auto setWithoutV = recursive();
        restoreElements(justV, connectedEdges);

        return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
    } else {
        return runIndependentSetDegree2();
    }
}

std::vector<NetworKit::node> Mis5IndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }

    auto v = getMaximumDegreeNode();
    if (graph->degree(v) >= 3) {
        auto neighborsPlus = getNeighborsPlus(v);
        EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
        removeElements(neighborsPlus);
        auto setWithV = recursive();
        setWithV.push_back(v);
        restoreElements(neighborsPlus, connectedEdges);

        auto justV = std::vector<NetworKit::node>{v};
        connectedEdges = getConnectedEdges(justV);
        removeElements(justV);
        auto setWithoutV = recursive();
        restoreElements(justV, connectedEdges);
        return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
    } else {
        return runIndependentSetDegree2();
    }
}

std::vector<NetworKit::node> MeasureAndConquerIndependentSet::recursive() {
    if (graph->isEmpty()) {
        return {};
    }
    std::vector<bool> visited(graph->upperNodeIdBound(), false);
    Koala::Traversal::DFSFrom(
        *graph, *graph->nodeRange().begin(),
        [&](auto v) { visited[v] = true; }, [&](auto v) { return true; });
    std::vector<NetworKit::node> component, theRest, allVertices;
    graph->forNodes([&](NetworKit::node u) {
        if (visited[u]) {
            component.push_back(u);
        } else {
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

    for (auto &nodeV : nodeNeighbors) {
        for (auto &nodeW : nodeNeighbors) {
            if (nodeV.index != nodeW.index) {
                if (std::includes(
                        nodeV.neighbors.begin(), nodeV.neighbors.end(),
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
    }  // vertices degree 1 don't exist at this point

    // A node v is foldable if N(v) contains no anti-triangles.
    // If none of the conditions above holds, and (3) there is a foldable node v of degree
    // d(v) <= 3, or a foldable node of degree d(v) = 4 with at most three anti-edges in
    // N(v), the algorithm selects one such node v of minimum degree and folds it
    for (auto v : graph->nodeRange()) {
        if (graph->degree(v) <= 4) {
            std::vector<NetworKit::node> vNeighbors = getNeighbors(v);
            EdgeSet inducedEdges = getInducedEdges(vNeighbors);
            bool shallFold = false;
            switch (graph->degree(v)) {
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
                // this also works: shallFold = true;
                if (inducedEdges.empty()) {
                    auto u1 = vNeighbors[0];
                    auto u2 = vNeighbors[1];
                    auto u12 = u1;
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
                } else {
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
                    if (inducedEdges.size() == 1) {  // other cases removed by domination
                        shallFold = true;
                    } else {
                        std::cout << "CAN'T BE HERE !!! (b)" << std::endl;
                    }
                }
                break;
            }
            case 4: {
                if (6 - inducedEdges.size() <= 3) {  // at most 3 anti-edges from K4
                    break;
                }
                if (inducedEdges.size() >= 4) {  // no triangle possible with 4 edges
                    shallFold = true;
                } else if (inducedEdges.size() == 3) {
                    // FOLD <=> no induced vertex has degree 3
                    bool hasDegree3 = false;
                    std::map<NetworKit::node, int> counter;
                    for (auto e : inducedEdges) {
                        ++counter[e.u], ++counter[e.v];
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
                    NetworKit::node oldU1, oldU2, u12;
                };
                std::vector<FoldingNode> foldingNodes;  // create one for each antiedge
                std::vector<std::set<NetworKit::node>> vNeighborsNeighbors(
                    graph->upperNodeIdBound());
                std::set<NetworKit::node> vNeighborsPlusSet(vNeighbors.begin(), vNeighbors.end());
                vNeighborsPlusSet.insert(v);

                for (auto u : vNeighbors) {
                    std::vector<NetworKit::node> uNeighborsVect = getNeighbors(u);
                    std::set<NetworKit::node> uNeighborsSet(
                        uNeighborsVect.begin(), uNeighborsVect.end());
                    std::set_difference(
                        uNeighborsSet.begin(), uNeighborsSet.end(),
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

    auto v = getMaximumDegreeNode();
    auto neighborsPlus = getNeighborsPlus(v);
    EdgeSet connectedEdges = getConnectedEdges(neighborsPlus);
    removeElements(neighborsPlus);
    auto setWithV = recursive();
    setWithV.push_back(v);
    restoreElements(neighborsPlus, connectedEdges);

    auto vWithMirrors = getMirrors(v);
    vWithMirrors.push_back(v);
    connectedEdges = getConnectedEdges(vWithMirrors);
    removeElements(vWithMirrors);
    auto setWithoutV = recursive();
    restoreElements(vWithMirrors, connectedEdges);

    return setWithV.size() > setWithoutV.size() ? setWithV : setWithoutV;
}

}  // namespace Koala
