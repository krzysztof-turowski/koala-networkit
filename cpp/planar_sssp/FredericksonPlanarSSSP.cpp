#include "planar_sssp/FredericksonPlanarSSSP.hpp"
#include "planar_sssp/PlanarUtils.hpp"

#include <set>
#include <unordered_map>
#include <utility>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <stdexcept>
#include <cmath>

using graphDivision_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

    NetworKit::node rootNode = 0;

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
        std::vector<std::vector<NetworKit::node>> components = cc.getComponents();

        for (auto c : Cbis) {

            for (auto x : G.neighborRange(c)) {
                NetworKit::count componentNumber = cc.componentOfNode(x);
                if (components[componentNumber].back() != c) {
                    components[cc.componentOfNode(x)].push_back(c);
                }
            }
        }
        return components;
    }

    std::vector<std::vector<NetworKit::node>> createConnectedSets(NetworKit::Graph& subGraph) {
        auto separator = findSeparator(subGraph);

        return postProcesing(separator, subGraph);
    }

    std::vector<std::unordered_set<NetworKit::node>> makeSuitableGraphDivision(const NetworKit::Graph& Graph, int r) {
        std::vector<std::vector<NetworKit::node>> smallSets;
        std::vector<std::vector<NetworKit::node>> smallAndAloneSets;
        std::vector<std::vector<NetworKit::node>> result;

        std::queue<std::vector<NetworKit::node>> queue;
        queue.push(std::vector<NetworKit::node> { Graph.nodeRange().begin(), Graph.nodeRange().end() });

        while (!queue.empty()) {
            auto nodeSet = queue.front();
            queue.pop();

            std::unordered_set<NetworKit::node> set{ nodeSet.begin(), nodeSet.end() };
            auto subGraph = NetworKit::GraphTools::subgraphFromNodes(Graph, set);
            auto sets = createConnectedSets(subGraph);
            for (auto& s : sets) {
                if (s.size() > r) {
                    queue.emplace(s);
                }
                else {
                    smallSets.push_back(s);
                }
            }
        }

        std::vector<std::vector<NetworKit::count>> componentsPerNode(Graph.numberOfNodes());
        std::vector<int> used(Graph.numberOfNodes(), 0);

        for (int i = 0; i < smallSets.size(); i++) {
            if (smallSets[i].size() > r / 2) {
                result.push_back(smallSets[i]);
                continue;
            }
            for (auto& x : smallSets[i]) {
                componentsPerNode[x].push_back(i);
            }
        }

        for (auto components : componentsPerNode) {
            assert(components.size() > 0 && components.size() < 4);

            if (components.size() == 1) {
                auto c = components[0];
                if (used[c]) continue;
                if (smallSets[c].size() > r / 2) {
                    result.push_back(smallSets[c]);
                }
                else {
                    smallAndAloneSets.push_back(smallSets[c]);
                }
            };


            if (components.size() == 2) {
                auto c0 = components[0];
                auto c1 = components[1];

            }

            if (components.size() == 3) {

            }
        }

        //TODO finish section 3
    }


    //FIND_CLUSTERS from [F1]
    graphDivision_t findClusterResult;
    std::vector<int> visited;

    std::vector<NetworKit::node> csearch(NetworKit::node v, int z, NetworKit::Graph& Graph) {
        if (visited[v]) {
            return {};
        }
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

    graphDivision_t findClusters(NetworKit::Graph& Graph, int z) {
        findClusterResult.clear();
        visited.assign(Graph.numberOfNodes(), 0);

        auto csearchResult = csearch(*Graph.nodeRange().begin(), z, Graph);

        for (auto node : csearchResult) {
            findClusterResult.back().push_back(node);
        }

        return findClusterResult;
    }

    graphDivision_t getDivisionOfFromClusters(graphDivision_t& clusters, NetworKit::Graph& Graph, int r) {
        //build new shrinked graph based on clusters
        NetworKit::Graph ShrinkedGraph(clusters.size() + 1);

        Graph.forNodes([&](auto v) {
            for (auto& u : Graph.neighborRange(v)) {
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

        graphDivision_t result(division.size(), std::vector<NetworKit::node> {});

        //Un Shrink the graph
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

    void fixGraphDivision(graphDivision_t& division, NetworKit::Graph& Graph) {
        //TODO implement infering boundy nodes
    }

    //create suitable r-division quickly.
    graphDivision_t  findSuitableRDivision(NetworKit::Graph Graph, int r) {
        int rSqrt = sqrt(r);

        auto clusters = findClusters(Graph, rSqrt);
        std::vector<NetworKit::count> clusterNum(Graph.numberOfNodes());
        for (int i = 0; i < clusters.size(); i++) {
            for (auto node : clusters[i]) {
                clusterNum[node] = i;
            }
        }

        graphDivision_t division = getDivisionOfFromClusters(clusters, Graph, r);
        //infer boundry nodes - whatever that means...
        fixGraphDivision(division, Graph);

        //TODO
        //run make Suitable Graph Divison from given regions
    }

    void FredericksonPlanarSSSP::run() {

        printGraph(graph);

        normal_graph = convertToMaxDegree3(graph);
        printGraph(normal_graph);



        std::cout << "------------------------------------------" << std::endl;
        auto [a, b, c] = findSeparator(normal_graph);

        for (auto x : a) {
            std::cout << x << " ";
        }std::cout << std::endl;
        for (auto x : b) {
            std::cout << x << " ";
        }std::cout << std::endl;
        for (auto x : c) {
            std::cout << x << " ";
        }std::cout << std::endl;

        std::vector<std::unordered_set<NetworKit::node>> Sets;
        std::queue<std::unordered_set<NetworKit::node>> setQueue;


        hasRun = true;
        distanceToTarget = 5;
        return;
    }

} /* namespace Koala */
