#include "planar_sssp/FredericksonPlanarSSSP.hpp"
#include "planar_sssp/PlanarUtils.hpp"
#include "planar_sssp/TopologyBasedHeap.hpp"

#include <set>
#include <unordered_map>
#include <utility>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <stdexcept>
#include <cmath>
#include <algorithm>

using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

    NetworKit::node rootNode = 0;
    std::vector<int> isBoundry;
    const int c = 10;

    NetworKit::Graph convertToMaxDegree3(NetworKit::Graph& G) {
        planar_embedding_t embedding = findPlanarEmbeding(G);
        NetworKit::count result_size = G.numberOfNodes();

        for (auto [i, node_embeding] : embedding) {
            if (node_embeding.size() > 3) {
                result_size += node_embeding.size() - 1;
            }
        }

        if (result_size == G.numberOfNodes()) {
            return G;
        }

        std::vector<std::unordered_map<size_t, size_t>> vertexMaps(
            G.numberOfNodes(), std::unordered_map<size_t, size_t>());

        size_t emptyVertex = G.numberOfNodes();
        NetworKit::Graph result(result_size, true, false);

        for (size_t i = 0; i < embedding.size(); ++i) {
            if (embedding[i].size() > 3) {
                size_t t = embedding[i][0];
                vertexMaps[i][t] = i;
                vertexMaps[i][emptyVertex] = 0;

                result.addEdge(i, emptyVertex, 0);
                for (int j = 1; j < embedding[i].size(); ++j) {
                    t = embedding[i][j];
                    vertexMaps[i][t] = emptyVertex;
                    if (j < embedding[i].size() - 1) {
                        result.addEdge(emptyVertex, emptyVertex + 1, 0);
                    }
                    emptyVertex++;
                }

                result.addEdge(emptyVertex - 1, i, 0);
            }
        }

        for (const auto& edge : G.edgeWeightRange()) {
            NetworKit::node u = edge.u;
            NetworKit::node v = edge.v;
            NetworKit::edgeweight ew = edge.weight;

            size_t s = v;
            size_t t = u;
            if (vertexMaps[u].find(v) != vertexMaps[u].end()) {
                t = vertexMaps[u][v];
            }
            if (vertexMaps[v].find(u) != vertexMaps[v].end()) {
                s = vertexMaps[v][u];
            }
            result.addEdge(s, t, ew);
        }

        return result;
    }


    struct Separator {
        std::unordered_set<NetworKit::node> A;
        std::unordered_set<NetworKit::node> B;
        std::unordered_set<NetworKit::node> C;
    };

    void printGraph(const NetworKit::Graph& graph) {
        std::cout << "Graph structure" << std::endl;
        for (auto edge : graph.edgeWeightRange()) {
            std::cout << edge.u << " - " << edge.v << " weight: " << edge.weight << std::endl;
        }
    }

    void insideDFS(NetworKit::node u, planar_embedding_t& graph, std::vector<bool>& visited, std::unordered_set<NetworKit::node>& inside) {
        if (visited[u]) return;
        inside.insert(u);
        visited[u] = true;
        for (auto v : graph[u]) {
            insideDFS(v, graph, visited, inside);
        }
        return;
    }

    std::pair<std::unordered_set<NetworKit::node>, std::unordered_set<NetworKit::node>>
        count_nodes(NetworKit::Edge edge, planar_embedding_t& graph, std::vector<NetworKit::node>& parent) {
        //TODO: improve cycle selection to be linear
        std::unordered_set<NetworKit::node> inside, cycle;

        std::vector<bool> visited(graph.size());
        NetworKit::node last;
        NetworKit::node current = edge.u;
        NetworKit::node rootRightChild;

        while (current != -1) {
            cycle.insert(current);
            visited[current] = true;
            current = parent[current];
        }
        current = edge.v;
        while (parent[current] != rootNode) {
            cycle.insert(current);
            visited[current] = true;
            current = parent[current];
        }
        rootRightChild = current;
        cycle.insert(current);
        visited[current] = true;

        // LEFT
        current = edge.u;
        last = edge.v;
        while (current != rootNode) {
            int start = std::find(graph[current].begin(), graph[current].end(), last) - graph[current].begin();
            int end = std::find(graph[current].begin(), graph[current].end(), parent[current]) - graph[current].begin();
            for (int i = (start + 1) % graph[current].size(); i < end; i = (i + 1) % graph[current].size()) {
                insideDFS(graph[current][i], graph, visited, inside);
            }
            last = current;
            current = parent[current];
        }

        // ROOT
        int start = std::find(graph[rootNode].begin(), graph[rootNode].end(), last) - graph[rootNode].begin();
        int end = std::find(graph[rootNode].begin(), graph[rootNode].end(), rootRightChild) - graph[rootNode].begin();
        for (int i = (start + 1) % graph[rootNode].size(); i < end; i = (i + 1) % graph[rootNode].size()) {
            insideDFS(graph[rootNode][i], graph, visited, inside);
        }

        return { inside, cycle };
    }

    Separator findSeparator(NetworKit::Graph& G) {
        NetworKit::Graph maximalGraph = makeMaximalPlanar(G);
        NetworKit::count numOfNodes = maximalGraph.numberOfNodes();

        std::vector<NetworKit::node> parent(numOfNodes, -1);
        std::vector<NetworKit::count> depth(numOfNodes, -1);
        std::vector<bool> visited(numOfNodes, 0);
        std::vector<NetworKit::Edge> treeEdges;
        std::vector<NetworKit::Edge> nonTreeEdges;

        std::queue<std::pair<NetworKit::node, NetworKit::count>> queue;

        auto planarEmbeding = findPlanarEmbeding(maximalGraph);
        parent[rootNode] = -1;
        depth[rootNode] = 0;
        visited[rootNode] = true;
        queue.push({ rootNode,0 });

        while (!queue.empty()) {
            auto [node, d] = queue.front();
            queue.pop();

            for (auto n : planarEmbeding[node]) {
                if (visited[n]) {
                    if (depth[node] <= depth[n]) {
                        nonTreeEdges.push_back({ node, n });
                    }
                    continue;
                }
                parent[n] = node;
                visited[n] = true;
                depth[n] = d + 1;
                treeEdges.push_back({ node, n });
                queue.push({ n, d + 1 });
            }
        }

        std::unordered_set<NetworKit::node> inside, cycle, outside;
        for (auto [u, v] : nonTreeEdges) {
            std::tie(inside, cycle) = count_nodes({ u, v }, planarEmbeding, parent);
            NetworKit::count outsideCount = numOfNodes - inside.size() - cycle.size();

            assert(inside.size() + cycle.size() <= numOfNodes);

            if (outsideCount * 3 <= 2 * numOfNodes && inside.size() * 3 <= 2 * numOfNodes) {
                maximalGraph.forNodes([&](NetworKit::node n) {
                    if (inside.find(n) == inside.end() && cycle.find(n) == cycle.end()) {
                        outside.insert(n);
                    }
                    });

                return { inside, outside, cycle };
            }
        }

        throw std::runtime_error("Acording to Lipton and Tarjan Lemma 2 Should not have happened!!!");
    }


    std::vector<std::vector<NetworKit::node>> postProcesing(Separator& sep, NetworKit::Graph& G) {
        auto [A, B, C] = sep;
        std::unordered_set<NetworKit::node> Cprime;
        std::unordered_set<NetworKit::node> Cbis;
        NetworKit::Graph copyG(G);

        for (auto c : C) {
            bool isBis = true;
            for (auto x : G.neighborRange(c)) {
                if (A.find(x) != A.end() || B.find(x) != B.end()) {
                    isBis = false;
                    Cprime.insert(c);
                    break;
                }
            }

            if (isBis) {
                Cbis.insert(c);
                for (auto x : G.neighborRange(c)) {
                    copyG.removeEdge(c, x);
                }
            }
        }
        NetworKit::ConnectedComponents cc(copyG);
        cc.run();
        nodeSubsets_t components = cc.getComponents();

        for (auto c : Cbis) {
            isBoundry[c] = 1;
            for (auto x : G.neighborRange(c)) {
                NetworKit::count componentNumber = cc.componentOfNode(x);
                if (components[componentNumber].back() != c) {
                    components[cc.componentOfNode(x)].push_back(c);
                }
            }
        }
        return components;
    }

    nodeSubsets_t createConnectedSets(NetworKit::Graph& subGraph) {
        auto separator = findSeparator(subGraph);

        return postProcesing(separator, subGraph);
    }

    bool canComponentBeMerged(int componentSize, int boundryNodeNum, int r, int sqr) {
        if (componentSize < r / 2 && boundryNodeNum < c * sqr) {
            return true;
        }
        return false;
    }

    //It requires isBoundry to be set correctly
    nodeSubsets_t makeSuitableGraphDivisionFromQueue(std::queue<std::vector<NetworKit::node>>& queue, NetworKit::Graph& Graph, int r) {
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
            }

            std::unordered_set<NetworKit::node> set{ nodeSet.begin(), nodeSet.end() };
            auto subGraph = NetworKit::GraphTools::subgraphFromNodes(Graph, set);
            auto sets = createConnectedSets(subGraph);
            for (auto& s : sets) {
                numOfBoundryNodes = 0;
                for (auto node : s) {
                    if (isBoundry[node]) {
                        numOfBoundryNodes++;
                    }
                }
                if (s.size() > r || numOfBoundryNodes > c * sqr) {
                    queue.emplace(std::move(s));
                }
                else {
                    smallSets.push_back(std::move(s));
                }
            }
        }

        //merging to small components
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

            //importat to merge smaller component to the bigger
            if (smallSets[c1].size() < smallSets[c2].size()) std::swap(c1, c2);
            isSetUsed[c2] = true;
            int goodBoundryNodes = 0;
            for (auto node : smallSets[c2]) {
                auto& nodeComponents = componentsPerNode[node];
                if (std::find(nodeComponents.begin(), nodeComponents.end(), c1) == nodeComponents.end()) {
                    smallSets[c1].push_back(node);
                    if (nodeComponents.size() > 1) {
                        goodBoundryNodes++;
                    }
                    std::replace(nodeComponents.begin(), nodeComponents.end(), c2, c1);
                }
                else {//boundry between c1 and c2
                    goodBoundryNodes += nodeComponents.size() == 2 ? -1 : 1;
                    nodeComponents.erase(find(nodeComponents.begin(), nodeComponents.end(), c2));
                }

            }
            numOfBoundryNodes[c1] += goodBoundryNodes;
            assert(numOfBoundryNodes[c1] >= 0);
            return c1;
            };

        //mergin components with common boundry nodes
        for (auto& components : componentsPerNode) {
            assert(components.size() > 0 && components.size() < 4);
            if (components.size() == 1) continue; //interior node

            tryMergeComponents(components[0], components[1]);
            if (components.size() == 2) {// boundry node of two regions
                continue;
            }
            // boundy node of threee regions
            tryMergeComponents(components[0], components[2]);
            tryMergeComponents(components[1], components[2]);
        }

        std::vector<std::tuple<int, int, int>> regionsToMerge;

        for (int i = 0; i < smallSets.size(); i++) {
            if (isSetUsed[i]) continue;
            std::unordered_set<int> regionNeighbours;
            for (auto node : smallSets[i]) {
                for (auto c : componentsPerNode[node]) {
                    regionNeighbours.insert(c);
                }
            }
            if (regionNeighbours.size() > 2) {
                // regions with more than two neighbours can be moved straight to the result
                // we can delay that to simplify the logic as the sets are still good to use
                continue;
            }

            assert(regionNeighbours.size() > 0); //there should NOT be a region without any neighbours

            int first, second;
            auto it = regionNeighbours.begin();
            first = *it;
            it++;
            second = it == regionNeighbours.end() ? -1 : *it;

            regionsToMerge.push_back({ first, second, i });
        }

        sort(regionsToMerge.begin(), regionsToMerge.end());//TODO: swap std::sort for an linear or nsqrt(logn) sort

        //greedy merging of regions with the same sets of either one or two regions.

        int lastComponent = -1;
        std::pair<int, int> lastPair = { -1,-1 };
        for (auto [f, s, i] : regionsToMerge) {
            if (!isComponentMergable(i)) continue;
            if (lastComponent = -1) {
                lastComponent = i;
                lastPair = { f, s };
                continue;
            }
            if (std::make_pair(f, s) != lastPair) {
                lastComponent = i;
                lastPair = { f,s };
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
            lastPair = { f,s };
        }

        for (int i = 0; i < smallSets.size(); i++) {
            if (isSetUsed[i]) continue;
            result.push_back(std::move(smallSets[i]));
        }

        return result;
    }

    // we need two different starting points for algorithm extra abstract layer allows for that
    nodeSubsets_t makeSuitableGraphDivision(NetworKit::Graph& Graph, int r) {
        isBoundry.assign(Graph.numberOfNodes(), 0);
        std::queue<std::vector<NetworKit::node>> queue;
        auto nodeRange = Graph.nodeRange();
        queue.push(std::vector<NetworKit::node> { nodeRange.begin(), nodeRange.end() });

        return makeSuitableGraphDivisionFromQueue(queue, Graph, r);
    }


    //FIND_CLUSTERS from [F1]
    nodeSubsets_t findClusterResult;
    std::vector<int> visited;

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
        }
        else {
            findClusterResult.push_back(clust);
            return {};
        }
    }

    nodeSubsets_t findClusters(NetworKit::Graph& Graph, int z) {
        findClusterResult.clear();
        visited.assign(Graph.numberOfNodes(), 0);

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
        //build new shrinked graph based on clusters
        NetworKit::Graph ShrinkedGraph(clusters.size() + 1);

        Graph.forNodes([&](auto v) {
            for (auto u : Graph.neighborRange(v)) {
                if (clusterNum[v] != clusterNum[u]) {
                    ShrinkedGraph.addEdge(clusterNum[v], clusterNum[u]);
                }
            }
            });
        ShrinkedGraph.removeMultiEdges();

        auto division = makeSuitableGraphDivision(ShrinkedGraph, r);

        std::vector<int> regionOfNode(ShrinkedGraph.numberOfNodes(), -1);
        std::vector<int> countRegionsOfNode(ShrinkedGraph.numberOfNodes(), 0);

        for (int i = 0; i < division.size(); i++) {
            for (auto& node : division[i]) {
                countRegionsOfNode[node] += 1;
                regionOfNode[node] = i;
            }
        }

        nodeSubsets_t result(division.size(), std::vector<NetworKit::node> {});

        //Unshrink the graph
        for (int i = 0; i < countRegionsOfNode.size(); i++) {
            assert(regionOfNode[i] > -1);
            //boundry vertex
            if (countRegionsOfNode[i] >= 2) {
                result.push_back({});
                for (auto node : clusters[i]) {
                    result.back().push_back(node);
                }
            }
            else {//interior
                for (auto node : clusters[i]) {
                    result[regionOfNode[i]].push_back(node);
                }
            }
        }

        return result;
    }

    void fixGraphDivision(nodeSubsets_t& division, NetworKit::Graph& Graph) {
        // The procedure is not described in the paper "Infer boundary vertices and slightly expanded(sic!) regions
        // that share these vertices". So I guess we can do it in a greedy fashion.

        std::vector<std::array<int, 3> > region(Graph.numberOfNodes(), { -1,-1,-1 });
        for (int i = 0; i < division.size(); i++) {
            for (auto& node : division[i]) {
                region[node][0] = i;
            }
        }

        Graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
            if (isBoundry[u] && isBoundry[v]) return;
            if (isBoundry[u] || isBoundry[v]) {
                if (isBoundry[v]) std::swap(u, v);
                if (region[u][1] == -1) {
                    region[u][1] = region[v][0];
                }
                if (region[u][1] != region[v][0]) {
                    region[u][2] = region[v][0];
                }
                return;
            }
            if (region[u][0] == region[v][0]) return;

            isBoundry[u] = 1;
            region[u][1] = region[v][0];
            });
    }

    //create suitable r-division quickly.
    nodeSubsets_t  findSuitableRDivision(NetworKit::Graph& Graph, int r) {
        int rSqrt = sqrt(r);

        auto clusters = findClusters(Graph, rSqrt);

        nodeSubsets_t division = getDivisionFromClusters(clusters, Graph, r);

        isBoundry.assign(Graph.numberOfNodes(), 0);
        fixGraphDivision(division, Graph);

        std::queue<std::vector<NetworKit::node>> queue;
        for (auto& region : division) {
            queue.push(std::move(region));
        }

        return makeSuitableGraphDivisionFromQueue(queue, Graph, r);
    }

    using pairDistance_t = std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;;
    pairDistance_t  getDistancebetweenBoundryNodesLevel2(NetworKit::Graph Graph, nodeSubsets_t division) {
        std::vector<pairDistance_t> resultPerRegion;
        pairDistance_t result;

        std::vector<int> regionCount(Graph.numberOfNodes(), 0);

        for (int i = 0; i < division.size(); i++) {
            for (auto node : division[i]) {
                regionCount[node]++;
            }
        }
        for (int i = 0; i < Graph.numberOfNodes(); i++) {
            if (regionCount[i] > 1) {
                isBoundry[i] = 1;
            }
        }

        for (int r = 0; r < division.size(); r++) {
            auto& region = division[r];
            auto subGraph = NetworKit::GraphTools::subgraphFromNodes(Graph,
                std::unordered_set<NetworKit::node> {region.begin(), region.end()});
            std::vector<NetworKit::node> boundry;

            for (auto node : region) {
                if (isBoundry[node]) {
                    boundry.push_back(node);
                }
            }

            for (int i = 0; i < boundry.size(); i++) {
                NetworKit::MultiTargetDijkstra dij(subGraph, boundry[i], boundry.begin(), boundry.end());
                dij.run();
                auto map = dij.getTargetIndexMap();
                auto distances = dij.getDistances();
                for (int j = 0; j < boundry.size(); j++) {
                    resultPerRegion[r][{boundry[i], boundry[j]}] = distances[map[boundry[j]]];
                    resultPerRegion[r][{boundry[j], boundry[i]}] = distances[map[boundry[j]]];
                }
            }
        }

        for (auto pairs : resultPerRegion) {
            for (auto [p, dist] : pairs) {
                if (result.contains(p)) {
                    result[p] = std::min(result[p], dist);
                }
                else {
                    result[p] = dist;
                }
            }
        }
        return result;
    }


    void mopUp(std::vector<int>& shortestDistance, nodeSubsets_t& regions, int s, int t) {

    }

    int mainThrust(NetworKit::Graph& Graph, nodeSubsets_t& regions, pairDistance_t& distances, int s, int t) {
        std::vector<int> shortestDistance(Graph.numberOfNodes(), -1);

        TopologyHeap heap(Graph, regions, distances, s);

        while (!heap.empty()) {
            auto [distance, node] = heap.top();
            shortestDistance[node] = distance;
            heap.closeNode(node);
        }

        mopUp(shortestDistance, regions, s, t);

        //TODO improve the mopUp function to calculate distance to all the nodes
        return shortestDistance[t];
    }

    void FredericksonPlanarSSSP::run() {
        printGraph(graph);
        normal_graph = convertToMaxDegree3(graph);
        printGraph(normal_graph);

        isBoundry.assign(normal_graph.numberOfNodes(), 0);
        int r = log(normal_graph.numberOfNodes());

        auto division = findSuitableRDivision(normal_graph, r);
        isBoundry.assign(normal_graph.numberOfNodes(), 0);

        auto distances = getDistancebetweenBoundryNodesLevel2(normal_graph, division);

        distanceToTarget = mainThrust(normal_graph, division, distances, source, target);

        hasRun = true;
        distanceToTarget = 5;
        return;
    }

} /* namespace Koala */
