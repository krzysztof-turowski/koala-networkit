#include "planar_sssp/SuitableRDivision.hpp"

#include <algorithm>
#include <cmath>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "planar_sssp/PlanarUtils.hpp"

using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {
NetworKit::node rootNode = 0;
std::vector<NetworKit::node> parent;
std::vector<NetworKit::count> depth;
std::vector<bool> visited;
int c;

struct Separator {
    std::unordered_set<NetworKit::node> A;
    std::unordered_set<NetworKit::node> B;
    std::unordered_set<NetworKit::node> C;
};

void insideDFS(NetworKit::node u, planar_embedding_t& graph,
               std::unordered_set<NetworKit::node>& inside) {
    if (visited[u]) return;
    inside.insert(u);
    visited[u] = true;
    for (auto v : graph[u]) {
        insideDFS(v, graph, inside);
    }
    return;
}

std::pair<std::unordered_set<NetworKit::node>, std::unordered_set<NetworKit::node>> count_nodes(
    NetworKit::Edge edge, planar_embedding_t& graph) {
    // TODO: improve cycle selection to be linear
    std::unordered_set<NetworKit::node> inside, cycle;

    for (auto& p : graph) {
        visited[p.first] = 0;
    }

    NetworKit::node last;
    NetworKit::node current = edge.u;
    NetworKit::node rootRightChild;
    NetworKit::node localRootNode = -1;

    while (current != -1) {
        cycle.insert(current);
        visited[current] = true;
        current = parent[current];
    }
    current = edge.v;
    while (!cycle.contains(parent[current])) {
        cycle.insert(current);
        visited[current] = true;
        current = parent[current];
    }
    rootRightChild = current;
    cycle.insert(current);
    visited[current] = true;
    localRootNode = parent[current];
    current = parent[localRootNode];
    while (current != -1) {
        cycle.erase(current);
        visited[current] = false;
        current = parent[current];
    }

    // LEFT
    current = edge.u;
    last = edge.v;
    while (current != localRootNode) {
        int start =
            std::find(graph[current].begin(), graph[current].end(), last) - graph[current].begin();
        int end = std::find(graph[current].begin(), graph[current].end(), parent[current]) -
                  graph[current].begin();
        for (int i = (start + 1) % graph[current].size(); i != end;
             i = (i + 1) % graph[current].size()) {
            insideDFS(graph[current][i], graph, inside);
        }
        last = current;
        current = parent[current];
    }

    // ROOT
    int start = std::find(graph[localRootNode].begin(), graph[localRootNode].end(), last) -
                graph[localRootNode].begin();
    int end = std::find(graph[localRootNode].begin(), graph[localRootNode].end(), rootRightChild) -
              graph[localRootNode].begin();
    for (int i = (start + 1) % graph[localRootNode].size(); i != end;
         i = (i + 1) % graph[localRootNode].size()) {
        insideDFS(graph[localRootNode][i], graph, inside);
    }

    return {inside, cycle};
}

Separator findSeparator(NetworKit::Graph& G) {
    rootNode = *(G.nodeRange().begin());
    NetworKit::Graph maximalGraph = makeMaximalPlanar(G);
    NetworKit::count numOfNodes = maximalGraph.numberOfNodes();

    for (auto node : maximalGraph.nodeRange()) {
        parent[node] = -1;
        depth[node] = -1;
        visited[node] = 0;
    }
    std::vector<NetworKit::Edge> treeEdges;
    std::vector<NetworKit::Edge> nonTreeEdges;

    std::queue<std::pair<NetworKit::node, NetworKit::count>> queue;

    auto planarEmbeding = findPlanarEmbeding(maximalGraph);
    parent[rootNode] = -1;
    depth[rootNode] = 0;
    visited[rootNode] = true;
    queue.push({rootNode, 0});

    while (!queue.empty()) {
        auto [node, d] = queue.front();
        queue.pop();

        for (auto n : planarEmbeding[node]) {
            if (visited[n]) {
                if (depth[node] <= depth[n]) {
                    nonTreeEdges.push_back({node, n});
                }
                continue;
            }
            parent[n] = node;
            visited[n] = true;
            depth[n] = d + 1;
            treeEdges.push_back({node, n});
            queue.push({n, d + 1});
        }
    }

    std::unordered_set<NetworKit::node> inside, cycle, outside;
    for (auto [u, v] : nonTreeEdges) {
        std::tie(inside, cycle) = count_nodes({u, v}, planarEmbeding);
        NetworKit::count outsideCount = numOfNodes - inside.size() - cycle.size();

        assert(inside.size() + cycle.size() <= numOfNodes);

        if (outsideCount * 3 <= 2 * numOfNodes && inside.size() * 3 <= 2 * numOfNodes) {
            maximalGraph.forNodes([&](NetworKit::node n) {
                if (inside.find(n) == inside.end() && cycle.find(n) == cycle.end()) {
                    outside.insert(n);
                }
            });

            for (auto [a, b] : G.edgeRange()) {
                if (cycle.contains(a) || cycle.contains(b)) continue;

                // no edge that connect inside with outside
                assert((inside.contains(a) && inside.contains(b)) ||
                       (outside.contains(a) && outside.contains(b)));
            }
            assert(inside.size() + outside.size() + cycle.size() == G.numberOfNodes());

            return {inside, outside, cycle};
        }
    }

    throw std::runtime_error("Acording to Lipton and Tarjan Lemma 2 Should not have happened!!!");
}

std::vector<std::vector<NetworKit::node>> postProcesing(Separator& sep, NetworKit::Graph& G,
                                                        std::vector<int>& isBoundry) {
    auto [A, B, C] = sep;
    std::unordered_set<NetworKit::node> Cprime;
    std::unordered_set<NetworKit::node> Cbis;
    NetworKit::Graph copyG(G);

    for (auto c : C) {
        bool isPrime = true;
        for (auto x : G.neighborRange(c)) {
            if (A.find(x) != A.end() || B.find(x) != B.end()) {
                isPrime = false;
                Cbis.insert(c);
                copyG.removeNode(c);
                break;
            }
        }

        if (isPrime) {
            Cprime.insert(c);
        }
    }
    NetworKit::ConnectedComponents cc(copyG);
    cc.run();
    nodeSubsets_t components = cc.getComponents();

    // the procedure was vaguely described in the paper
    // what's important is that the process of moving vertices from Cbis to coneccted componets is
    // greedy it's easy to prove that bellow prcess will created actuall connected subsets it
    // requires observation about structure of subgraph Cbis. Cbis is a uninion of unconnected
    // paths.

    std::unordered_map<int, int> movedNodeMap;
    std::unordered_set<NetworKit::node> boundryNodes;
    for (auto c : Cbis) {
        std::unordered_set<int> neighbourComponents;
        for (auto x : G.neighborRange(c)) {
            if (Cbis.contains(x)) {
                if (movedNodeMap.find(x) != movedNodeMap.end()) {
                    neighbourComponents.insert(movedNodeMap[x]);
                }
                continue;
            }
            neighbourComponents.insert(cc.componentOfNode(x));
        }

        if (neighbourComponents.size() == 1) {
            movedNodeMap[c] = *(neighbourComponents.begin());
            components[movedNodeMap[c]].push_back(c);
        } else {
            boundryNodes.insert(c);
        }
    }

    std::unordered_map<int, int> newBoundryNode;
    for (auto node : boundryNodes) {
        isBoundry[node] = 1;
        newBoundryNode[node] = 1;
        for (auto x : G.neighborRange(node)) {
            int component;

            if (Cbis.contains(x)) {
                if (movedNodeMap.find(x) == movedNodeMap.end()) {
                    continue;
                }
                component = movedNodeMap[x];
            } else {
                component = cc.componentOfNode(x);
            }

            if (components[component].back() != node) {
                components[component].push_back(node);
            }
            movedNodeMap[node] = component;
        }
    }

    // Assert will be removed later

    std::unordered_map<int, std::vector<int>> regionsOfNode;
    for (int i = 0; i < components.size(); i++) {
        for (auto node : components[i]) {
            regionsOfNode[node].push_back(i);
        }
    }

    for (auto& [k, v] : regionsOfNode) {
        int countNeighbours = 0;
        for (auto nei : G.neighborRange(k)) {
            countNeighbours++;
        }
        assert(v.size() > 0);
        assert(v.size() <= countNeighbours);
        if (newBoundryNode[k]) assert(v.size() > 1);
    }

    for (auto [u, v] : G.edgeRange()) {
        bool hasCommon = std::any_of(regionsOfNode[u].begin(), regionsOfNode[u].end(), [&](int x) {
            return std::find(regionsOfNode[v].begin(), regionsOfNode[v].end(), x) !=
                   regionsOfNode[v].end();
        });
        assert(hasCommon && "edge must be in at least one region fully");
        // in other words baundry nodes are in all regions they are connected to
    }

    // Assert ends here

    return components;
}

nodeSubsets_t createConnectedSets(NetworKit::Graph& subGraph, std::vector<int>& isBoundry) {
    nodeSubsets_t result;
    NetworKit::ConnectedComponents cc(subGraph);
    cc.run();
    nodeSubsets_t components = cc.getComponents();
    if (components.size() == 1) {  // subgraph is a connected componet
        auto separator = findSeparator(subGraph);
        return postProcesing(separator, subGraph, isBoundry);
    }

    int largestComponentNumber = -1;
    for (int i = 0; i < components.size(); i++) {
        if (components[i].size() * 3 > subGraph.numberOfNodes() * 2) {
            largestComponentNumber = i;
        } else {
            result.push_back(std::move(components[i]));
        }
    }
    if (largestComponentNumber == -1) {  // all components ale smaller thatn 2/3 of number of nodes
                                         // no need for using separater theorem
        return result;
    }

    // finding separator of the largest component (it size exceeds 2/3 off all nodes)
    std::unordered_set<NetworKit::node> largestComponent(components[largestComponentNumber].begin(),
                                                         components[largestComponentNumber].end());
    NetworKit::Graph connectedGraph =
        NetworKit::GraphTools::subgraphFromNodes(subGraph, largestComponent);
    auto separator = findSeparator(connectedGraph);

    auto conectedSubsets = postProcesing(separator, connectedGraph, isBoundry);
    for (int i = 0; i < conectedSubsets.size(); i++) {
        result.push_back(std::move(conectedSubsets[i]));
    }

    return result;
}

bool canComponentBeMerged(int componentSize, int boundryNodeNum, int r, int sqr) {
    if (componentSize < r / 2 && boundryNodeNum < c * sqr) {
        return true;
    }
    return false;
}

void fixDivision(nodeSubsets_t& division, NetworKit::Graph& Graph) {
    std::vector<std::vector<int>> componentsOfNode(Graph.numberOfNodes());

    for (int i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            componentsOfNode[node].push_back(i);
        }
    }

    for (auto node : Graph.nodeRange()) {
        if (componentsOfNode[node].size() > 3) {
            int componentToRemove = -1;

            for (auto c : componentsOfNode[node]) {
                int countNeighbors = 0;
                int TheNei = -1;
                for (auto nei : Graph.neighborRange(node)) {
                    if (std::find(componentsOfNode[nei].begin(), componentsOfNode[nei].end(), c) !=
                        componentsOfNode[nei].end()) {
                        countNeighbors++;
                        TheNei = nei;
                    }
                }

                if (countNeighbors == 0) {
                    componentToRemove = c;
                    break;
                }
                if (countNeighbors == 1) {
                    std::vector<int> restComponents;
                    for (auto cc : componentsOfNode[node]) {
                        if (c != cc) {
                            restComponents.push_back(cc);
                        }
                    }
                    bool hasCommon =
                        std::any_of(restComponents.begin(), restComponents.end(), [&](int x) {
                            return std::find(componentsOfNode[TheNei].begin(),
                                             componentsOfNode[TheNei].end(),
                                             x) != componentsOfNode[TheNei].end();
                        });
                    if (hasCommon) {
                        componentToRemove = c;
                        break;
                    }
                }
            }

            division[componentToRemove].erase(std::find(division[componentToRemove].begin(),
                                                        division[componentToRemove].end(), node));
        }
    }
}

nodeSubsets_t makeSuitableGraphDivisionFromQueue(std::queue<std::vector<NetworKit::node>>& queue,
                                                 NetworKit::Graph& Graph, int r,
                                                 std::vector<int>& isBoundry) {
    int sqr = sqrt(r);
    nodeSubsets_t smallSets;
    nodeSubsets_t result;

    // dividing graph into small overlaping sets of nodes
    while (!queue.empty()) {
        auto nodeSet = std::move(queue.front());
        queue.pop();

        int numOfBoundryNodes = 0;
        for (auto node : nodeSet) {
            if (isBoundry[node]) {
                numOfBoundryNodes++;
            }
        }
        if (nodeSet.size() < r && numOfBoundryNodes < c * sqr) {
            smallSets.push_back(std::move(nodeSet));
            continue;
        }

        std::unordered_set<NetworKit::node> nodesForSubGraph{nodeSet.begin(), nodeSet.end()};
        auto subGraph = NetworKit::GraphTools::subgraphFromNodes(Graph, nodesForSubGraph);
        auto sets = createConnectedSets(subGraph, isBoundry);

        for (auto& s : sets) {
            numOfBoundryNodes = 0;
            for (auto node : s) {
                if (isBoundry[node]) {
                    numOfBoundryNodes++;
                }
            }
            if (s.size() >= r || numOfBoundryNodes >= c * sqr) {
                queue.emplace(std::move(s));
            } else {
                smallSets.push_back(std::move(s));
            }
        }
    }

    fixDivision(smallSets, Graph);

    assert_division(smallSets, Graph);

    // merging to small components
    std::vector<std::vector<NetworKit::count>> componentsPerNode(Graph.numberOfNodes());
    std::vector<int> isSetUsed(smallSets.size(), 0);
    std::vector<int> numOfBoundryNodes(smallSets.size(), 0);
    std::vector<int> used(Graph.numberOfNodes(), 0);

    for (int i = 0; i < smallSets.size(); i++) {
        for (auto& x : smallSets[i]) {
            componentsPerNode[x].push_back(i);
        }
    }
    for (int i = 0; i < componentsPerNode.size(); i++) {
        if (componentsPerNode[i].size() > 1) {
            for (auto c : componentsPerNode[i]) {
                numOfBoundryNodes[c] += 1;
            }
        }
    }

    auto isComponentMergable = [&](int c) {
        if (isSetUsed[c]) return false;
        if (2 * smallSets[c].size() > r || 2 * numOfBoundryNodes[c] > c * sqr) return false;
        return true;
    };

    auto tryMergeComponents = [&](int c1, int c2) {
        if (!isComponentMergable(c1)) return -1;
        if (!isComponentMergable(c2)) return -1;

        // importat to merge smaller component to the bigger
        if (smallSets[c1].size() < smallSets[c2].size()) std::swap(c1, c2);
        isSetUsed[c2] = true;
        int goodBoundryNodes = 0;
        for (auto node : smallSets[c2]) {
            auto& nodeComponents = componentsPerNode[node];
            if (std::find(nodeComponents.begin(), nodeComponents.end(), c1) ==
                nodeComponents.end()) {
                smallSets[c1].push_back(node);
                if (nodeComponents.size() > 1) {
                    goodBoundryNodes++;
                }
                std::replace(nodeComponents.begin(), nodeComponents.end(), c2, c1);
            } else {  // boundry between c1 and c2
                goodBoundryNodes += nodeComponents.size() == 2 ? -1 : 1;
                nodeComponents.erase(find(nodeComponents.begin(), nodeComponents.end(), c2));
            }
        }
        smallSets[c2].clear();
        numOfBoundryNodes[c1] += goodBoundryNodes;
        assert(numOfBoundryNodes[c1] >= 0);
        return c1;
    };

    // mergin components with common boundry nodes
    for (auto& components : componentsPerNode) {
        if (components.size() == 1) continue;  // interior node

        int c0, c1, c2;
        if (components.size() == 3) {  // boundy node of threee regions
            int c0 = components[0], c1 = components[1], c2 = components[2];
            tryMergeComponents(c0, c1);
            tryMergeComponents(c0, c2);
            tryMergeComponents(c1, c2);
            continue;
        }
        // boundry node of two regions
        tryMergeComponents(components[0], components[1]);
    }

    std::vector<std::vector<NetworKit::node>> regionsToAssert;
    for (auto set : smallSets) {
        if (set.size() > 0) {
            regionsToAssert.push_back(set);
        }
    }
    assert_division(regionsToAssert, Graph);

    int countBoundryNodes = 0;
    for (auto node : Graph.nodeRange()) {
        assert(componentsPerNode[node].size() > 0);
        if (componentsPerNode[node].size() > 1) {
            countBoundryNodes++;
        }
    }

    std::vector<std::tuple<int, int, int>> regionsToMerge;

    for (int i = 0; i < smallSets.size(); i++) {
        if (smallSets[i].size() == 0) continue;
        if (!isComponentMergable(i)) continue;
        std::unordered_set<int> regionNeighbours;
        for (auto node : smallSets[i]) {
            for (auto c : componentsPerNode[node]) {
                if (c != i) {
                    regionNeighbours.insert(c);
                }
            }
        }
        if (regionNeighbours.size() > 2) {
            // regions with more than two neighbours can be moved straight to the result
            // we can delay that to simplify the logic as the sets are still good to use
            continue;
        }

        int first, second;
        auto it = regionNeighbours.begin();
        first = *it;
        it++;
        second = it == regionNeighbours.end() ? -1 : *it;

        regionsToMerge.push_back({first, second, i});
    }

    sort(regionsToMerge.begin(),
         regionsToMerge.end());  // TODO: swap std::sort for an linear or nsqrt(logn) sort

    // greedy merging of regions with the same sets of either one or two regions.

    int lastComponent = -1;
    std::pair<int, int> lastPair = {-1, -1};
    for (auto [f, s, i] : regionsToMerge) {
        if (!isComponentMergable(i)) continue;
        if (lastComponent = -1) {
            lastComponent = i;
            lastPair = {f, s};
            continue;
        }
        if (std::make_pair(f, s) != lastPair) {
            lastComponent = i;
            lastPair = {f, s};
            continue;
        }

        int merged = tryMergeComponents(lastComponent, i);
        if (!isComponentMergable(merged)) {
            lastComponent = -1;
            continue;
        }
        if (merged == lastComponent) {
            continue;
        }
        lastComponent = merged;
        lastPair = {f, s};
    }

    for (int i = 0; i < smallSets.size(); i++) {
        if (isSetUsed[i]) continue;
        result.push_back(std::move(smallSets[i]));
    }

    assert_division(result, Graph);
    return result;
}

// we need two different starting points for algorithm
nodeSubsets_t makeSuitableGraphDivision(NetworKit::Graph& Graph, int r) {
    std::vector<int> isBoundry(Graph.numberOfNodes(), 0);
    std::queue<std::vector<NetworKit::node>> queue;
    auto nodeRange = Graph.nodeRange();
    queue.push(std::vector<NetworKit::node>{nodeRange.begin(), nodeRange.end()});

    return makeSuitableGraphDivisionFromQueue(queue, Graph, r, isBoundry);
}

// FIND_CLUSTERS from [F1]
nodeSubsets_t findClusterResult;

std::vector<NetworKit::node> csearch(NetworKit::node v, int z, NetworKit::Graph& Graph) {
    if (visited[v]) {
        return {};
    }
    visited[v] = 1;
    std::vector<NetworKit::node> clust(1, v);

    for (auto neigh : Graph.neighborRange(v)) {
        auto neighClust = csearch(neigh, z, Graph);
        for (auto n : neighClust) {
            clust.push_back(n);
        }
    }
    if (clust.size() < z) {
        return clust;
    } else {
        findClusterResult.push_back(clust);
        return {};
    }
}

nodeSubsets_t findClusters(NetworKit::Graph& Graph, int z) {
    findClusterResult.clear();
    for (auto node : Graph.nodeRange()) {
        visited[node] = 0;
    }

    auto csearchResult = csearch(*Graph.nodeRange().begin(), z, Graph);

    for (auto node : csearchResult) {
        findClusterResult.back().push_back(node);
    }

    return findClusterResult;
}

nodeSubsets_t getDivisionFromClusters(nodeSubsets_t& clusters, NetworKit::Graph& Graph, int r) {
    std::vector<NetworKit::count> clusterNum(Graph.numberOfNodes());
    for (int i = 0; i < clusters.size(); i++) {
        for (auto node : clusters[i]) {
            clusterNum[node] = i;
        }
    }
    // build new shrinked graph based on clusters
    NetworKit::Graph ShrinkedGraph(clusters.size());

    Graph.forNodes([&](auto v) {
        for (auto u : Graph.neighborRange(v)) {
            if (clusterNum[v] != clusterNum[u]) {
                ShrinkedGraph.addEdge(clusterNum[v], clusterNum[u]);
            }
        }
    });
    ShrinkedGraph.removeMultiEdges();

    int sqr = sqrt(r);
    nodeSubsets_t division;

    std::vector<int> isBoundry(ShrinkedGraph.numberOfNodes(), 0);
    std::queue<std::vector<NetworKit::node>> queue;
    auto nodeRange = ShrinkedGraph.nodeRange();
    queue.push(std::vector<NetworKit::node>{nodeRange.begin(), nodeRange.end()});

    // dividing graph into small overlaping sets of nodes
    while (!queue.empty()) {
        auto nodeSet = std::move(queue.front());
        queue.pop();

        int numOfBoundryNodes = 0;
        for (auto node : nodeSet) {
            if (isBoundry[node]) {
                numOfBoundryNodes++;
            }
        }
        if (nodeSet.size() < r && numOfBoundryNodes < c * sqr) {
            division.push_back(std::move(nodeSet));
            continue;
        }

        std::unordered_set<NetworKit::node> nodesForSubGraph{nodeSet.begin(), nodeSet.end()};
        auto subGraph = NetworKit::GraphTools::subgraphFromNodes(ShrinkedGraph, nodesForSubGraph);
        auto sets = createConnectedSets(subGraph, isBoundry);

        for (auto& s : sets) {
            numOfBoundryNodes = 0;
            for (auto node : s) {
                if (isBoundry[node]) {
                    numOfBoundryNodes++;
                }
            }
            if (s.size() >= r || numOfBoundryNodes >= c * sqr) {
                queue.emplace(std::move(s));
            } else {
                division.push_back(std::move(s));
            }
        }
    }

    std::vector<int> regionOfNode(ShrinkedGraph.numberOfNodes(), -1);
    std::vector<int> countRegionsOfNode(ShrinkedGraph.numberOfNodes(), 0);

    for (int i = 0; i < division.size(); i++) {
        for (auto& node : division[i]) {
            countRegionsOfNode[node] += 1;
            regionOfNode[node] = i;
        }
    }

    nodeSubsets_t resultWithZeros(division.size(), std::vector<NetworKit::node>{});

    // Unshrink the graph
    for (int i = 0; i < countRegionsOfNode.size(); i++) {
        assert(regionOfNode[i] > -1);
        // boundry vertex
        if (countRegionsOfNode[i] >= 2) {
            resultWithZeros.push_back({});
            for (auto node : clusters[i]) {
                resultWithZeros.back().push_back(node);
            }
        } else {  // interior
            for (auto node : clusters[i]) {
                resultWithZeros[regionOfNode[i]].push_back(node);
            }
        }
    }

    nodeSubsets_t result;

    for (int i = 0; i < resultWithZeros.size(); i++) {
        if (resultWithZeros[i].size() > 0) {
            result.push_back(std::move(resultWithZeros[i]));
        }
    }

    return result;
}

nodeSubsets_t fixGraphDivision(nodeSubsets_t& division, NetworKit::Graph& Graph,
                               std::vector<int>& isBoundry) {
    // The procedure is not described in the paper "Infer boundary vertices and slightly
    // expanded(sic!) regions that share these vertices". So I guess we can do it in a greedy
    // fashion.

    std::vector<std::array<int, 3>> region(Graph.numberOfNodes(), {-1, -1, -1});
    for (int i = 0; i < division.size(); i++) {
        for (auto& node : division[i]) {
            region[node][0] = i;
        }
    }

    Graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        if (isBoundry[u] && isBoundry[v]) {
            return;
        }
        if (isBoundry[u] || isBoundry[v]) {
            if (isBoundry[v]) std::swap(u, v);
            for (auto r : region[u]) {
                if (r == region[v][0]) return;
            }
            region[u][2] = region[v][0];
            return;
        }
        if (region[u][0] == region[v][0]) return;

        isBoundry[u] = 1;
        region[u][1] = region[v][0];
    });

    Graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        if (isBoundry[u] && isBoundry[v]) {
            bool hasCommon = false;

            for (auto cu : region[u]) {
                if (cu == -1) break;
                for (auto cv : region[v]) {
                    if (cv == cu) hasCommon = true;
                }
            }
            if (!hasCommon) {
                if (region[u][2] == -1) {
                    region[u][2] = region[v][0];
                } else {
                    region[v][2] = region[u][0];
                }
            }
        }
    });

    nodeSubsets_t result(division.size());

    for (int i = 0; i < region.size(); i++) {
        for (auto reg : region[i]) {
            if (reg == -1) continue;
            result[reg].push_back(i);
        }
    }

    return result;
}

// create suitable r-division quickly.
nodeSubsets_t findSuitableRDivision(NetworKit::Graph& Graph, int r, int constC) {
    // values used globally in dfs
    parent.resize(Graph.numberOfNodes());
    depth.resize(Graph.numberOfNodes());
    visited.resize(Graph.numberOfNodes());
    c = constC;
    int rSqrt = sqrt(r);

    auto clusters = findClusters(Graph, rSqrt);

    nodeSubsets_t division = getDivisionFromClusters(clusters, Graph, r);

    for (auto& region : division) {
        assert(region.size() > 0);
    }

    std::vector<int> isBoundry(Graph.numberOfNodes(), 0);
    division = fixGraphDivision(division, Graph, isBoundry);

    assert_division(division, Graph);

    std::queue<std::vector<NetworKit::node>> queue;
    for (auto& region : division) {
        queue.push(std::move(region));
    }
    return makeSuitableGraphDivisionFromQueue(queue, Graph, r, isBoundry);
}
}  // namespace Koala
