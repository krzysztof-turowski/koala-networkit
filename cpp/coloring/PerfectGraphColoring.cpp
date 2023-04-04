/*
 * PerfectGraphColoring.cpp
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

#include <coloring/PerfectGraphColoring.hpp>
#include <graph/GraphTools.hpp>

#include "perfect/commons.h"

// Input: Graph of n vertices and m edges. Arrays from and to of size m, from[i]->to[i] is i'th edge.
// Nodes start at 1.
double get_theta(int n, int m, int *from, int *to) {
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
    constraints = static_cast<struct constraintmatrix*>(malloc((m + 2) * sizeof(struct constraintmatrix)));
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
        constraints[i].blocks = static_cast<struct sparseblock*>(malloc(sizeof(struct sparseblock)));
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
    return getTheta(Koala::GraphTools::toComplement(graph), get_vector(graph));
}

std::vector<int> get_maximum_stable_set(const NetworKit::Graph &graph) {
    std::vector<int> keep_nodes(get_vector(graph));
    int theta = getTheta(graph, keep_nodes);
    graph.forNodes([&](NetworKit::node v) {
        keep_nodes[v] = 0;
        if (getTheta(graph, keep_nodes) != theta) {
            keep_nodes[v] = 1;
        }
    });
    return keep_nodes;
}

std::vector<int> get_maximum_clique(const NetworKit::Graph &graph) {
  return get_maximum_stable_set(Koala::GraphTools::toComplement(graph));
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

std::vector<NetworKit::node> PerfectGraphColoring::get_stable_set_intersecting_all_maximum_cliques() {
    auto K2 = get_maximum_clique(*graph);
    int omega2 = std::accumulate(K2.begin(), K2.end(), 0);
    while (true) {
        std::vector<int> S2 = get_stable_set_intersecting_maximum_cliques_2(K2);
        NetworKit::Graph subgraph2 = *graph;
        for (auto v : S2) {
            subgraph2.removeNode(v);
        }
        int omega4 = getOmega(subgraph2);
        if (omega4 < omega2) {
            std::vector<NetworKit::node> out;
            for (auto i : S2) {
                out.push_back(i);
            }
            return out;
        }
        auto clique3 = get_maximum_clique(subgraph2);
        for (int i = 0; i < clique3.size(); i++) {
            K2[i] += clique3[i];
        }
    }
}

std::vector<int> PerfectGraphColoring::get_stable_set_intersecting_maximum_cliques_2(const std::vector<int> &K) {
    printf("K2: ");
    for(auto v : K) printf("%d ", v);
    printf("\n");
    std::vector<int> c(K.size());
    std::partial_sum(K.begin(), K.end(), c.begin());
    NetworKit::Graph auxiliary_graph(c.back(), false, false);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        for (int ni = u == 0 ? 0 : c[u - 1]; ni < c[u]; ni++) {
            for (int nj = v == 0 ? 0 : c[v - 1]; nj < c[v]; nj++) {
                auxiliary_graph.addEdge(ni, nj);
            }
        }
    });
    auto stable_set = get_maximum_stable_set(auxiliary_graph);
    // printf("stable_set: ");
    // for(auto v : stable_set) printf("%d ", v);
    // printf("\n");
    std::vector<int> out;
    for (int i = 0, index = 0; i < stable_set.size(); i++) {
        if (stable_set[i]) {
            while (index < c.size() && c[index] <= i) {
                index++;
            }
            if (out.empty() || out.back() != index) {
                out.push_back(index);
            }
        }
    }
    // printf("stable_set_out: ");
    // for(auto v : out) printf("%d ", v);
    // printf("\n");
    return out;
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
