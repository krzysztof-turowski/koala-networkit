/*
 * MinimumSpanningTree.cpp
 *
 *  Created on: 07.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <mst/MinimumSpanningTree.hpp>

#include <random>
#include <ranges>
#include <queue>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/graph/GraphTools.hpp>

#include "mst/MSTV.hpp"
#include "mst/KktMST.hpp"
  
std::default_random_engine generator;
std::uniform_int_distribution<int> distribution(0, 1);

namespace Koala {

MinimumSpanningTree::MinimumSpanningTree(
        NetworKit::Graph &graph) : graph(std::make_optional(graph)) {
    tree = std::make_optional(NetworKit::GraphTools::copyNodes(graph));
}

const NetworKit::Graph& MinimumSpanningTree::getForest() const {
    assureFinished();
    return *tree;
}

void KruskalMinimumSpanningTree::run() {
    hasRun = true;
    std::vector<NetworKit::WeightedEdge> sorted_edges(
        graph->edgeWeightRange().begin(), graph->edgeWeightRange().end());
    Aux::Parallel::sort(sorted_edges.begin(), sorted_edges.end());
    NetworKit::UnionFind union_find(graph->upperNodeIdBound());
    for (const auto &e : sorted_edges) {
        if (union_find.find(e.u) != union_find.find(e.v)) {
            tree->addEdge(e.u, e.v, e.weight);
            union_find.merge(e.u, e.v);
        }
    }
}

void PrimMinimumSpanningTree::run() {
    hasRun = true;
    std::priority_queue<std::pair<NetworKit::edgeweight, NetworKit::node>> queue;
    queue.push(std::make_pair(0, *(graph->nodeRange().begin())));
    std::unordered_map<NetworKit::node, NetworKit::WeightedEdge> previous;
    while (!queue.empty()) {
        auto v = queue.top().second;
        queue.pop();
        if (!tree->isIsolated(v)) {
            continue;
        }
        const auto &e = previous.find(v);
        if (e != previous.end()) {
            tree->addEdge(e->second.u, e->second.v, e->second.weight);
        }
        graph->forNeighborsOf(v, [&](NetworKit::node u, NetworKit::edgeweight weight) {
            if (tree->isIsolated(u)) {
                queue.push(std::make_pair(-weight, u));
                if (!previous.count(u) || previous[u].weight > weight) {
                    previous[u] = NetworKit::WeightedEdge(u, v, weight);
                }
            }
        });
    }
}

void BoruvkaMinimumSpanningTree::run() {
    hasRun = true;
    NetworKit::UnionFind union_find(graph->upperNodeIdBound());
    iterate(*graph, *tree, union_find, std::numeric_limits<NetworKit::count>::max());
}

void BoruvkaMinimumSpanningTree::iterate(
        NetworKit::Graph &G, NetworKit::Graph &F, NetworKit::UnionFind &union_find,
        NetworKit::count steps) {
    using Edge = std::pair<NetworKit::node, NetworKit::node>;
    std::map<Edge, Edge> E;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({{std::min(u, v), std::max(u, v)}, {u, v}});
    });
    while (G.numberOfNodes() > 1 && G.numberOfEdges() > 0 && steps-- > 0) {
        G.forNodes([&](NetworKit::node x) {
            auto e = *std::min_element(
                G.weightNeighborRange(x).begin(), G.weightNeighborRange(x).end(),
                [&](const auto &e1, const auto &e2) { return e1.second < e2.second; });
            NetworKit::node u_prim = union_find.find(x), v_prim = union_find.find(e.first);
            if (u_prim == v_prim) {
                return;
            }
            if (x > e.first) {
                std::swap(x, e.first);
            }
            union_find.merge(u_prim, v_prim);
            const auto &[u, v] = E[{x, e.first}];
            F.addEdge(u, v, e.second);
        });
        std::map<NetworKit::node, std::vector<NetworKit::WeightedEdge>> first_pass;
        for (auto e : G.edgeWeightRange()) {
            const NetworKit::node &u_prim = union_find.find(e.u), &v_prim = union_find.find(e.v);
            if (u_prim < v_prim) {
                std::swap(e.u, e.v);
                first_pass[v_prim].emplace_back(std::move(e));
            } else if (u_prim > v_prim) {
                first_pass[u_prim].emplace_back(std::move(e));
            }
        }
        std::map<NetworKit::node, std::vector<NetworKit::WeightedEdge>> second_pass;
        for (const auto &[v_prim, E_v] : first_pass | std::views::reverse) {
            for (const auto &e : E_v) {
                second_pass[union_find.find(e.u)].emplace_back(std::move(e));
            }
        }
        G.removeAllEdges();
        for (const auto &v : G.nodeRange()) {
            if (union_find.find(v) != v) {
                G.removeNode(v);
            }
        }
        for (const auto &[u_prim, E_u] : second_pass) {
            for (auto left = E_u.begin(); left != E_u.end();) {
                const auto &v_prim = union_find.find(left->v);
                auto right = std::find_if(
                    left + 1, E_u.end(),
                    [&](const auto &e) { return union_find.find(e.v) != v_prim; });
                const auto &e = *std::min_element(
                    left, right,
                    [&](const auto &e1, const auto &e2) { return e1.weight < e2.weight; });
                G.addEdge(u_prim, v_prim, e.weight);
                E[{u_prim, v_prim}] = E[{std::min(e.u, e.v), std::max(e.u, e.v)}];
                left = right;
            }
        }
    }
}

void KargerKleinTarjanMinimumSpanningTree::run() {
    hasRun = true;
    recurse(*graph, *tree);
}

void KargerKleinTarjanMinimumSpanningTree::recurse(NetworKit::Graph &G, NetworKit::Graph &F) {
    NetworKit::UnionFind union_find(G.upperNodeIdBound());
    while (true) {
        iterate(G, F, union_find, 2);
        if (G.numberOfEdges() == 0) {
            return;
        }
        auto subgraph(NetworKit::GraphTools::copyNodes(G));
        discard_random_edges(G, subgraph);
        auto subforest(NetworKit::GraphTools::copyNodes(subgraph));
        recurse(subgraph, subforest);
        remove_heavy_edges(G, subforest);
    }
}

void KargerKleinTarjanMinimumSpanningTree::discard_random_edges(
        NetworKit::Graph &G, NetworKit::Graph &subgraph) {
    G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (distribution(generator)) {
            subgraph.addEdge(u, v, w);
        }
    });
    auto connected_components = NetworKit::ConnectedComponents(subgraph);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (auto i = 1; i < connected_components.numberOfComponents(); i++) {
        subgraph.addEdge(
            components[0][0], components[i][0], std::numeric_limits<NetworKit::edgeweight>::max());
    }
}

void KargerKleinTarjanMinimumSpanningTree::remove_heavy_edges(
        NetworKit::Graph &G, NetworKit::Graph &subforest) {
    MST::removeF_HeavyEdges(G, subforest);
}

} /* namespace Koala */
