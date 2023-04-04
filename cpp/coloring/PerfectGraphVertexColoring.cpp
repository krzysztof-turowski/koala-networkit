/*
 * PerfectGraphVertexColoring.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <map>
#include <optional>
#include <tuple>

extern "C" {
#include <declarations.h>
}

#include <coloring/PerfectGraphVertexColoring.hpp>
#include <graph/GraphTools.hpp>

constexpr int CACHE_LIMIT = 6;

// Graph of n vertices and m edges. Arrays from and to of size m, from[i]->to[i] is i'th edge.
// Nodes start at 1.
double compute_theta(int n, int m, int *from, int *to) {
    double pobj, dobj, *y, *a;
    struct blockmatrix C, X, Z;
    struct constraintmatrix *constraints;

    C.nblocks = 1;
    C.blocks = static_cast<blockrec*>(malloc(2 * sizeof(struct blockrec)));
    C.blocks[1].blockcategory = MATRIX;
    C.blocks[1].blocksize = n;
    C.blocks[1].data.mat = static_cast<double*>(malloc(n * n * sizeof(double)));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            C.blocks[1].data.mat[ijtok(i, j, n)] = 1.0;
        }
    }
    a = static_cast<double*>(malloc((m + 2)* sizeof(double)));
    a[1] = 1.0;
    constraints = static_cast<struct constraintmatrix*>(
        malloc((m + 2) * sizeof(struct constraintmatrix)));
    constraints[1].blocks = static_cast<struct sparseblock*>(malloc(sizeof(struct sparseblock)));
    constraints[1].blocks->blocknum = constraints[1].blocks->constraintnum = 1;
    constraints[1].blocks->numentries = constraints[1].blocks->blocksize = n;
    constraints[1].blocks->next = constraints[1].blocks->nextbyblock = NULL;
    constraints[1].blocks->entries = static_cast<double*>(malloc((n + 1) * sizeof(double)));
    constraints[1].blocks->iindices = static_cast<int*>(malloc((n + 1) * sizeof(int)));
    constraints[1].blocks->jindices = static_cast<int*>(malloc((n + 1) * sizeof(int)));
    for (int i = 1; i <= n; i++) {
        constraints[1].blocks->entries[i] = 1.0;
        constraints[1].blocks->iindices[i] = constraints[1].blocks->jindices[i] = i;
    }
    for (int i = 2; i <= m + 1; i++) {
        a[i] = 0.0;
        constraints[i].blocks = static_cast<struct sparseblock*>(
            malloc(sizeof(struct sparseblock)));
        constraints[i].blocks->blocknum = constraints[i].blocks->numentries = 1;
        constraints[i].blocks->blocksize = n;
        constraints[i].blocks->constraintnum = i;
        constraints[i].blocks->next = constraints[i].blocks->nextbyblock = NULL;
        constraints[i].blocks->entries = static_cast<double*>(malloc(2 * sizeof(double)));
        constraints[i].blocks->iindices = static_cast<int*>(malloc(2 * sizeof(int)));
        constraints[i].blocks->jindices = static_cast<int*>(malloc(2 * sizeof(int)));
        constraints[i].blocks->entries[1] = 1.0;
        int start = from[i - 2], finish = to[i - 2];
        if (start > finish) {
            std::swap(start, finish);
        }
        constraints[i].blocks->iindices[1] = start;
        constraints[i].blocks->jindices[1] = finish;
    }
    initsoln(n, m + 1, C, a, constraints, &X, &y, &Z);
    easy_sdp(n, m + 1, C, a, constraints, 0.0, &X, &y, &Z, &pobj, &dobj);
    free_prob(n, m + 1, C, a, constraints, X, y, Z);
    return (dobj + pobj) / 2;
}

std::vector<int> get_vector(const NetworKit::Graph &graph) {
    std::vector<int> out(graph.upperNodeIdBound(), 0);
    graph.forNodes([&](NetworKit::node v) {
        out[v] = 1;
    });
    return out;
}

namespace Koala {

void PerfectGraphVertexColoring::run() {
    omega = get_omega(*graph);
    for (int color = 1; graph->numberOfNodes() > 0; color++) {
        for (const auto &v : get_stable_set_intersecting_all_maximum_cliques()) {
            colors[v] = color;
            graph->removeNode(v);
        }
    }
    hasRun = true;
}

void PerfectGraphVertexColoring::check() const {
    assureFinished();
    int chi = std::max_element(
        std::begin(colors), std::end(colors),
        [] (const auto &a, const auto &b) { return a.second < b.second; })->second;
    assert(omega == chi);
}

std::vector<int> PerfectGraphVertexColoring::get_stable_set_intersecting_all_maximum_cliques() {
    auto cliques = get_maximum_clique(*graph);
    int omega = std::accumulate(cliques.begin(), cliques.end(), 0);
    while (true) {
        auto stable_set = get_maximum_weighted_stable_set(*graph, cliques);
        NetworKit::Graph subgraph = *graph;
        for (auto v : stable_set) {
            subgraph.removeNode(v);
        }
        if (get_omega(subgraph) < omega) {
            return stable_set;
        }
        auto clique = get_maximum_clique(subgraph);
        for (int i = 0; i < clique.size(); i++) {
            cliques[i] += clique[i];
        }
    }
}

int PerfectGraphVertexColoring::get_theta(
        const NetworKit::Graph &graph, const std::vector<int> &keep_nodes) {
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
    static std::map<decltype(e), int> CACHE;
    if (indices.back() <= CACHE_LIMIT && CACHE.count(e) > 0) {
        return CACHE[e];
    }
    double theta = compute_theta(indices.back(), from.size(), from.data(), to.data());
    int theta_int = theta + 0.5;
    double eps = 0.3;
    if (abs(theta - theta_int) > eps) {
        throw std::logic_error("Non-integer theta for a perfect graph: " + std::to_string(theta));
    }
    if (indices.back() <= CACHE_LIMIT) {
        CACHE[e] = theta_int;
    }
    return theta_int;
}

int PerfectGraphVertexColoring::get_omega(const NetworKit::Graph &graph) {
    return get_theta(Koala::GraphTools::toComplement(graph), get_vector(graph));
}

std::vector<int> PerfectGraphVertexColoring::get_maximum_clique(const NetworKit::Graph &graph) {
    return get_maximum_stable_set(Koala::GraphTools::toComplement(graph));
}

std::vector<int> PerfectGraphVertexColoring::get_maximum_stable_set(
        const NetworKit::Graph &graph) {
    std::vector<int> keep_nodes(get_vector(graph));
    int theta = get_theta(graph, keep_nodes);
    graph.forNodes([&](NetworKit::node v) {
        keep_nodes[v] = 0;
        if (get_theta(graph, keep_nodes) != theta) {
            keep_nodes[v] = 1;
        }
    });
    return keep_nodes;
}

std::vector<int> PerfectGraphVertexColoring::get_maximum_weighted_stable_set(
        const NetworKit::Graph &graph, const std::vector<int> &weights) {
    std::vector<int> count(weights.size());
    std::partial_sum(weights.begin(), weights.end(), count.begin());
    NetworKit::Graph auxiliary_graph(count.back(), false, false);
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        for (int ni = u == 0 ? 0 : count[u - 1]; ni < count[u]; ni++) {
            for (int nj = v == 0 ? 0 : count[v - 1]; nj < count[v]; nj++) {
                auxiliary_graph.addEdge(ni, nj);
            }
        }
    });
    auto stable_set = get_maximum_stable_set(auxiliary_graph);
    std::vector<int> out;
    for (int i = 0, index = 0; i < stable_set.size(); i++) {
        if (stable_set[i]) {
            index = std::find_if(
                count.begin() + index, count.end(), [&](int v) { return v > i; }) - count.begin();
            if (out.empty() || out.back() != index) {
                out.push_back(index);
            }
        }
    }
    return out;
}

} /* namespace Koala */
