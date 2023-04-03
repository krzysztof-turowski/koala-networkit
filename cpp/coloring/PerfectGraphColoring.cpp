/*
 * PerfectGraphColoring.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <map>
#include <optional>
#include <tuple>

#include <coloring/PerfectGraphColoring.hpp>
#include <graph/GraphTools.hpp>

#include "perfect/commons.h"
#include "perfect/theta.h"

Graph generate_graph(std::optional<NetworKit::Graph> graph) {
    std::map<NetworKit::node, int> vertices;
    int nodes = 0;
    graph->forNodes([&](NetworKit::node v) {
        vertices[v] = nodes++;
    });
    std::vector<std::vector<int>> M(nodes);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        M[vertices[u]].push_back(vertices[v]), M[vertices[v]].push_back(vertices[u]);
    });
    return Graph(M);
}

NetworKit::Graph retrieve_graph(const Graph &G) {
    NetworKit::Graph graph(G.n, false, false);
    for (int i = 0; i < G.n; i++) {
        for (int j = i + 1; j < G.n; j++) {
          if (G.areNeighbours(i, j)) {
              graph.addEdge(i, j);
          }
        }
    }
    return graph;
}

int getTheta(const NetworKit::Graph &graph, const std::vector<int> &keep_nodes) {
    std::vector<int> indices(keep_nodes.size());
    std::partial_sum(keep_nodes.begin(), keep_nodes.end(), indices.begin());
    if (indices.empty() || indices.back() == 0) {
        return 0;
    }
    if (indices.back() == 1) {
        return 1;
    }

    std::vector<int> from, to;
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        if (keep_nodes[u] && keep_nodes[v]) {
            from.push_back(indices[u]), to.push_back(indices[v]);
        }
    });
    if (indices.back() == 2) {
        return 1 + (from.size() == 0);
    }

    auto e = std::make_tuple(indices.back(), from.size(), from, to);
    static std::map<std::tuple<int, int, std::vector<int>, std::vector<int>>, int> CACHE;
    if (CACHE.count(e) > 0) {
        return CACHE[e];
    }

    double theta = get_theta(indices.back(), from.size(), from.data(), to.data());
    int theta_int = theta + 0.5;
    double eps = 0.3;
    if (abs(theta - theta_int) > eps) {
        throw std::logic_error("Theta returned non-integer for a Perfect Graph: " + std::to_string(theta));
    }
    return CACHE[e] = theta_int;
}

std::vector<int> get_vector(const NetworKit::Graph &graph) {
    std::vector<int> out(graph.upperNodeIdBound(), 0);
    graph.forNodes([&](NetworKit::node v) {
        out[v] = 1;
    });
    return out;
}

int getOmega(const NetworKit::Graph &graph) {
    auto graph_complement = Koala::GraphTools::toComplement(graph);
    assert(graph.numberOfNodes() == graph_complement.numberOfNodes());
    return getTheta(graph_complement, get_vector(graph_complement));
}

int getOmega(const Graph &G) {
    NetworKit::Graph graph(retrieve_graph(G));
    auto graph_complement = Koala::GraphTools::toComplement(graph);
    return getTheta(graph_complement, get_vector(graph_complement));
}

std::vector<int> get_maximum_stable_set(const NetworKit::Graph &graph) {
    std::vector<int> keep_nodes(get_vector(graph));
    int theta = getTheta(graph, keep_nodes);
    std::vector<int> res;
    graph.forNodes([&](NetworKit::node v) {
        keep_nodes[v] = 0;
        if (getTheta(graph, keep_nodes) != theta) {
            keep_nodes[v] = 1;
            res.push_back(v);
        }
    });
    return res;
}

std::vector<int> get_maximum_clique(const Graph &G) {
  NetworKit::Graph graph(retrieve_graph(G));
  auto graph_complement = Koala::GraphTools::toComplement(graph);
  return get_maximum_stable_set(graph_complement);
}

namespace Koala {

int PerfectGraphColoring::get_omega() {
    return getOmega(*graph);
}

int PerfectGraphColoring::get_chi() {
    run();
    int out = 0;
    for (auto &[_, v] : colors) {
        out = std::max(out, v);
    }
    return out;
}

std::set<NetworKit::node> PerfectGraphColoring::get_stable_set_intersecting_all_maximum_cliques() {
    Graph G(generate_graph(graph));
    std::vector<int> vertices;
    graph->forNodes([&](NetworKit::node v) {
        vertices.push_back(v);
    });
    auto maximum_clique = get_maximum_clique(G);
    int omega = maximum_clique.size();
    std::vector<int> K(graph->numberOfNodes(), 0);
    for (auto v : maximum_clique) {
        K[v] = 1;
    }
    while (true) {
        std::vector<int> S = get_stable_set_intersecting_maximum_cliques(K);
        std::vector<int> compS = getComplementNodesVec(G.n, S);
        Graph Gprim = G.getInducedStrong(compS);
        if (getOmega(Gprim) < omega) {
            std::set<NetworKit::node> out;
            for (auto i : S) {
                out.insert(vertices[i]);
            }
            return out;
        } else {
            auto clique = get_maximum_clique(Gprim);
            for (auto v : clique) {
                K[compS[v]]++;
            }
        }
    }
}

std::vector<int> PerfectGraphColoring::get_stable_set_intersecting_maximum_cliques(const std::vector<int> &K) {
    std::vector<int> c(K.size());
    std::map<NetworKit::node, int> V;
    int nodes = 0;
    graph->forNodes([&](NetworKit::node v) {
        V[v] = nodes++;
    });
    std::partial_sum(K.begin(), K.end(), c.begin());
    NetworKit::Graph auxiliary_graph(c.back(), false, false);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        for (int ni = V[u] == 0 ? 0 : c[V[u] - 1]; ni < c[V[u]]; ni++) {
            for (int nj = V[v] == 0 ? 0 : c[V[v] - 1]; nj < c[V[v]]; nj++) {
                auxiliary_graph.addEdge(ni, nj);
            }
        }
    });
    std::vector<int> nSS = get_maximum_stable_set(auxiliary_graph);
    std::vector<int> ret;
    int wsk = 0;
    for (int nSSnode : nSS) {
        while (wsk < graph->numberOfNodes() && c[wsk] <= nSSnode) wsk++;
        if (ret.empty() || ret.back() != wsk) {
            ret.push_back(wsk);
        }
    }
    return ret;
}

void PerfectGraphColoring::run() {
    for (int color = 1; graph->numberOfNodes() > 0; color++) {
        for (const auto &v : get_stable_set_intersecting_all_maximum_cliques()) {
            colors[v] = color;
            graph->removeNode(v);
        }
    }
    hasRun = true;
}

} /* namespace Koala */
