#include "planar_sssp/FredericksonPlanarSSSP.hpp"
#include "planar_sssp/PlanarUtils.hpp"

#include <set>
#include <unordered_map>
#include <utility>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/BFS.hpp>
#include <stdexcept>

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
            bool canBeAdded = true;
            auto lastComponent = cc.componentOfNode(*(G.neighborRange(c).begin()));

            for (auto x : G.neighborRange(c)) {
                if (lastComponent != cc.componentOfNode(x)) {
                    canBeAdded = false;
                    break;
                }
                lastComponent = cc.componentOfNode(x);
            }

            if (canBeAdded) {
                Cbis.erase(c);
                components[cc.componentOfNode(*(G.neighborRange(c).begin()))].push_back(c);
            }
        }
        components.push_back(std::vector<NetworKit::node> {Cbis.begin(), Cbis.end()});

        return components;
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

    std::vector<std::vector<NetworKit::node>> createConnectedSets(NetworKit::Graph& subGraph) {
        auto separator = findSeparator(subGraph);

        return postProcesing(separator, subGraph);
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

    void FredericksonPlanarSSSP::run() {

        // NetworKit::Graph subgraph = NetworKit::GraphTools::subgraphFromNodes(graph, { 2,3 });

        // std::cout << subgraph.numberOfNodes() << std::endl;

        // for (auto x : subgraph.nodeRange()) {
        //     std::cout << x << " ";
        // }std::cout << std::endl;

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
