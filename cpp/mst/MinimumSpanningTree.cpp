/*
 * MinimumSpanningTree.cpp
 *
 *  Created on: 07.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <mst/MinimumSpanningTree.hpp>

#include <random>
#include <ranges>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <structures/Heap.hpp>

#include "mst/MSTV.hpp"
#include "mst/BoruvkaMST.hpp"
#include "mst/KktMST.hpp"

std::random_device device;
std::default_random_engine generator{device()};
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
    Heap<std::pair<NetworKit::edgeweight, NetworKit::node>> queue;
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
    NetworKit::Graph G(*graph);
    iterate(G, *tree, union_find, std::numeric_limits<NetworKit::count>::max(), true);
}

void initialize_branching_tree(
        NetworKit::Graph &B, std::unordered_map<NetworKit::node, NetworKit::node> &V_B,
        const NetworKit::Graph &G) {
    B = NetworKit::Graph(G.upperNodeIdBound(), true, true);
    for (const auto &v : G.nodeRange()) {
        V_B[v] = v;
    }
}

std::optional<NetworKit::Graph> BoruvkaMinimumSpanningTree::iterate(
        NetworKit::Graph &G, NetworKit::Graph &F, NetworKit::UnionFind &union_find,
        NetworKit::count steps, bool get_branching_tree) {
    using Edge = std::pair<NetworKit::node, NetworKit::node>;
    std::map<Edge, Edge> E;
    G.forEdges([&](const NetworKit::node u, const NetworKit::node v) {
        E.insert({std::minmax(u, v), {u, v}});
    });

    NetworKit::Graph B;
    std::unordered_map<NetworKit::node, NetworKit::node> V_B;
    if (get_branching_tree) {
        std::cout << "INITIALIZE" << std::endl;
        initialize_branching_tree(B, V_B, G);
    }
    while (G.numberOfNodes() > 1 && G.numberOfEdges() > 0 && steps-- > 0) {
        std::unordered_map<NetworKit::node, NetworKit::edgeweight> B_edges;
        G.forNodes([&](NetworKit::node x) {
            auto [y, w] = *std::min_element(
                G.weightNeighborRange(x).begin(), G.weightNeighborRange(x).end(),
                [](const auto &e1, const auto &e2) { return e1.second < e2.second; });
            if (get_branching_tree) {
                B_edges[x] = w;
            }
            NetworKit::node u_prim = union_find.find(x), v_prim = union_find.find(y);
            if (u_prim == v_prim) {
                return;
            }
            union_find.merge(u_prim, v_prim);
            const auto &[u, v] = E[std::minmax(x, y)];
            F.addEdge(u, v, w);
        });
        std::map<NetworKit::node, std::vector<NetworKit::WeightedEdge>> first_pass;
        for (auto e : G.edgeWeightRange()) {
            const NetworKit::node &u_prim = union_find.find(e.u), &v_prim = union_find.find(e.v);
            if (u_prim < v_prim) {
                first_pass[v_prim].emplace_back(std::move(e));
            } else if (u_prim > v_prim) {
                std::swap(e.u, e.v);
                first_pass[u_prim].emplace_back(std::move(e));
            }
        }
        std::map<NetworKit::node, std::vector<NetworKit::WeightedEdge>> second_pass;
        for (const auto &[v_prim, E_v] : (first_pass | std::views::reverse)) {
            for (const auto &e : E_v) {
                second_pass[union_find.find(e.u)].emplace_back(std::move(e));
            }
        }
        auto G_prim = NetworKit::GraphTools::copyNodes(G);
        for (const auto &v : G_prim.nodeRange()) {
            if (union_find.find(v) != v) {
                G_prim.removeNode(v);
            }
        }
        if (get_branching_tree) {
            std::unordered_map<NetworKit::node, NetworKit::node> V_B_next;
            for (const auto &v : G_prim.nodeRange()) {
                V_B_next[v] = B.addNode();
            }
            std::cout << "X-ADDED " << G_prim.numberOfNodes() << " " << B_edges.size() << std::endl;
            for (const auto &[v, w] : B_edges) {
                // std::cout << "EDGE " << V_B_next[union_find.find(v)] << " " << V_B[v] << " " << w << std::endl;
                B.addEdge(V_B_next[union_find.find(v)], V_B[v], w);
            }
            std::swap(V_B, V_B_next);
        }
        for (const auto &[u_prim, E_u] : second_pass) {
            auto left = E_u.begin();
            while (left != E_u.end()) {
                const auto &v_prim = union_find.find(left->v);
                auto right = std::find_if(
                    left + 1, E_u.end(),
                    [&](const auto &e) { return union_find.find(e.v) != v_prim; });
                const auto &e = *std::min_element(
                    left, right,
                    [](const auto &e1, const auto &e2) { return e1.weight < e2.weight; });
                assert(G_prim.hasNode(u_prim) && G_prim.hasNode(v_prim));
                G_prim.addEdge(u_prim, v_prim, e.weight);
                E[{u_prim, v_prim}] = E[std::minmax(e.u, e.v)];
                left = right;
            }
        }
        std::swap(G, G_prim);
    }
    return get_branching_tree ? std::make_optional(B) : std::nullopt;
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
        std::cout << std::endl << G.numberOfEdges() << " " << G.numberOfEdges() << std::endl;
        /*auto subgraph(NetworKit::GraphTools::copyNodes(G));
        discard_random_edges(G, union_find, subgraph);
        auto subforest(NetworKit::GraphTools::copyNodes(subgraph));
        recurse(subgraph, subforest);
        remove_heavy_edges(G, subforest);*/
        auto S = NetworKit::GraphTools::getContinuousNodeIds(G);
        auto G_prim = NetworKit::GraphTools::getCompactedGraph(G, S);
        std::cout << G_prim.numberOfEdges() << " " << G_prim.numberOfEdges() << std::endl;
        auto algorithm = MST::KktMST(G_prim);
        auto S_invert = NetworKit::GraphTools::invertContinuousNodeIds(S, G);
        for (auto [u, v] : algorithm.getSpanningTree().edgeRange()) {
            F.addEdge(u, v, G_prim.weight(u, v));
        }
        break;
    }
}

void KargerKleinTarjanMinimumSpanningTree::discard_random_edges(
        NetworKit::Graph &G, NetworKit::UnionFind &union_find, NetworKit::Graph &subgraph) {
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
    auto S = NetworKit::GraphTools::getContinuousNodeIds(G);
    auto G_prim = NetworKit::GraphTools::getCompactedGraph(G, S);
    auto subforest_prim = NetworKit::GraphTools::getCompactedGraph(subforest, S);
    auto S_invert = NetworKit::GraphTools::invertContinuousNodeIds(S, G);
    MST::removeF_HeavyEdges(G_prim, subforest_prim, G, subforest, S_invert);
}

class AugmentedGraph {
    NetworKit::Graph graph;
 public:
    AugmentedGraph(NetworKit::Graph &G) : graph(G) {
    }

    std::optional<NetworKit::node> getParent(NetworKit::node v) {
        if (graph.inNeighborRange(v).begin() == graph.inNeighborRange(v).end()) {
            return std::nullopt;
        }
        return *(graph.inNeighborRange(v).begin());
    }

    NetworKit::node getRoot() {
        return graph.upperNodeIdBound() - 1;
    }

    NetworKit::Graph& getTree() {
        return graph;
    }
};

void MinimumSpanningTree::check() const {
    assureFinished();
    NetworKit::Graph tree_copy(*tree);
    NetworKit::UnionFind union_find(graph->upperNodeIdBound());
    NetworKit::Graph branching_tree_base = *BoruvkaMinimumSpanningTree::iterate(
        tree_copy, tree_copy, union_find, std::numeric_limits<NetworKit::count>::max(), true);

    std::cout << "CHECK" << std::endl;
    // Boruvka runs in linear time on trees.
    MST::BoruvkaMST fbtOnMst(*tree, true, true);
    std::vector<std::tuple<NetworKit::node, NetworKit::node, NetworKit::edgeweight>> G_minus_M;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
         if (tree->hasEdge(u, v)) {
            G_minus_M.push_back({u, v, w});
         }
    });
    std::cout << "CHECK 2" << std::endl;
    const NetworKit::Graph& mst = *tree, &G = *graph;
    auto branching_tree2 = AugmentedGraph(branching_tree_base);
    auto lca2 = Koala::LCA(branching_tree2);
    auto &branching_tree = fbtOnMst.fullBranchingTree.value();
    std::map<NetworKit::node, NetworKit::node> B;
    std::cout << branching_tree_base.numberOfNodes() << " " << branching_tree_base.numberOfEdges() << std::endl;
    std::cout << branching_tree.getTree().numberOfNodes() << " " << branching_tree.getTree().numberOfEdges() << std::endl;
    for (auto v : graph->nodeRange()) {
        B[v] = v;
    }
    for (auto v : branching_tree.getTree().nodeRange()) {
        if (v == branching_tree.getTree().upperNodeIdBound() - 1) {
            break;
        }
        auto v_parent = branching_tree.getParent(v).value();
        auto u = B[v], u_parent = branching_tree2.getParent(u).value();
        B[v_parent] = u_parent;
        if (v < graph->numberOfNodes()) {
            assert(u == v);
        }
    }
    branching_tree_base.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (v == branching_tree_base.upperNodeIdBound() - 1) {
            assert(branching_tree_base.degreeIn(v) == 0);
        } else {
            assert(branching_tree_base.degreeIn(v) == 1);
        }
    });
    auto lca = Koala::LCA(branching_tree);
    std::vector<NetworKit::node> upper, lower;
    std::vector<NetworKit::node> upper2, lower2;
    for (auto [u, v, w] : G_minus_M) {
        auto lowestCommonAncestor = lca.query(u, v);
        upper.push_back(lowestCommonAncestor), upper.push_back(lowestCommonAncestor);
        lower.push_back(u), lower.push_back(v);
        
        auto lowestCommonAncestor2 = lca2.query(u, v);
        upper2.push_back(lowestCommonAncestor2), upper2.push_back(lowestCommonAncestor2);
        lower2.push_back(u), lower2.push_back(v);
    }

    std::cout << "CHECK 3" << std::endl;
    auto answers = MST::treePathMaxima(branching_tree, upper, lower);
    auto answers2 = MST::treePathMaxima(branching_tree2, upper2, lower2);
    for (NetworKit::index i = 0; i < answers.size(); i++) {
        assert(MST::edgeWeightToParent(answers2[i], branching_tree2) <= std::get<2>(G_minus_M[i / 2]));
        assert(MST::edgeWeightToParent(answers[i], branching_tree) <= std::get<2>(G_minus_M[i / 2]));
    }
    std::cout << "CHECK 4" << std::endl;
}

} /* namespace Koala */
