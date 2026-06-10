/*
 * MinimumSpanningTree.cpp
 *
 *  Created on: 07.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <ranges>

#include <algorithm>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>

#include "graph/GraphTools.hpp"
#include "mst/MinimumSpanningTree.hpp"
#include "structures/Heap.hpp"
#include "structures/LCA.hpp"

std::random_device device;
std::default_random_engine generator{device()};
std::uniform_int_distribution<int> distribution(0, 1);

class BranchingTree {
 private:
    NetworKit::Graph B;
    std::unordered_map<NetworKit::node, NetworKit::node> V_B;
    std::unordered_map<NetworKit::node, NetworKit::edgeweight> B_edges;

 public:
    void initialize(const NetworKit::Graph &G) {
        B = NetworKit::Graph(G.upperNodeIdBound(), true, true);
        for (const auto &v : G.nodeRange()) {
            V_B[v] = v;
        }
    }

    inline void addEdge(NetworKit::node v, NetworKit::edgeweight w) {
        B_edges[v] = w;
    }

    void update(const NetworKit::Graph &G, const NetworKit::UnionFind &union_find) {
        std::unordered_map<NetworKit::node, NetworKit::node> V_B_next;
        for (const auto &v : G.nodeRange()) {
            V_B_next[v] = B.addNode();
        }
        for (const auto &[v, w] : B_edges) {
            B.addEdge(V_B_next[union_find.find(v)], V_B[v], w);
        }
        std::swap(V_B, V_B_next);
        B_edges.clear();
    }

    NetworKit::Graph getTree() {
        return B;
    }
};

class AugmentedGraph {
    using h_set = uint64_t;  // bitset of depths. i-th bit corresponds to i-th depth.

    NetworKit::Graph graph;

    NetworKit::count n, height;
    // median: h_set -> element of h_set. Example: 0b1011 denotes {0,1,3} so median[0b1011] = 1
    std::vector<NetworKit::index> median;
    // depth: node -> depth.
    std::vector<NetworKit::count> depth;
    // D[u]: a set of depths of endpoints above u of query paths that contain u. D[root] = 0.
    std::vector<h_set> D;
    // P: depth -> node. DFS stack, used in `visit`.
    std::vector<NetworKit::node> P;
    // `L` and `L_next` combined form a linked list of "indices of queries" (query_ids).
    std::vector<NetworKit::index> L, L_next;

 public:
    explicit AugmentedGraph(NetworKit::Graph &G) : graph(G), n(G.numberOfNodes()) { }

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

    NetworKit::edgeweight getWeight(NetworKit::node u) {
        if (u == getRoot()) {
            return 0;
        }
        auto w = getTree().weight(getParent(u).value(), u);
        assert(w > 0);
        return w;
    }

    std::vector<NetworKit::node> getTreePathMaxima(
            const std::vector<NetworKit::node>& lower, const std::vector<NetworKit::node>& upper) {
        height = 0, depth.resize(n, 0), D.resize(n, 0);
        L.resize(n, NetworKit::none), L_next.resize(upper.size(), NetworKit::none);
        for (std::size_t i = 0; i < lower.size(); i++) {  // distribute queries to lower nodes.
            // L[u] - beginning of a linked list query_ids that have lower node `u`.
            // L[u] == index of the first query. L_next[L[u]] - index of the next query.
            // L[] is indexed by nodes and L_next[] is indexed by query_ids.
            // L[u] points to the beginning of the list for `u`. This list is stored in `L_next[]`.
            L_next[i] = L[lower[i]], L[lower[i]] = i;
        }
        initialize(getRoot(), 0, upper);
        P.resize(height + 1);
        median_table(height);
        std::vector<NetworKit::node> answer(upper.size(), NetworKit::none);
        visit(getRoot(), 0, upper, answer);
        return answer;
    }

 private:
    void initialize(
            NetworKit::node u, NetworKit::count u_depth,
            const std::vector<NetworKit::node> &upper) {
        depth[u] = u_depth;
        if (u_depth > height) {
            height = u_depth;
        }
        for (auto i = L[u]; i != NetworKit::none; i = L_next[i]) {
            D[u] |= 1 << depth[upper[i]];  // this is executed only for leaves of `fbt`.
        }
        getTree().forNeighborsOf(u, [&](NetworKit::node child) {
            initialize(child, u_depth + 1, upper);
            // exclude `u`. no-op for leaves and works recursively up to the root.
            D[u] |= (D[child] & ~(1 << u_depth));
        });
    }

    void median_table(NetworKit::count h) {
        // Fills a table of size 2^(h+1) whose entry in position i, for
        // i = 0, ..., 2^(h-1) - 1, is the median of the set represented by i.
        std::vector<h_set> T((1 << h) + 1);
        median.resize(1 << (h + 1));
        for (NetworKit::count s = 0; s <= h; s++) {
            for (NetworKit::count k = 0; k <= s; k++) {
                auto p = subsets(T, h - s, k, 0);
                auto q = subsets(T, s, k + 1, subsets(T, s, k, p));
                for (NetworKit::count i = 0; i < p; i++) {
                    auto b = (1 << (s + 1)) * T[i] + (1 << s);  // fixed high bits
                    for (auto j = p; j < q; j++) {
                        median[b + T[j]] = s;  // variable low bits
                    }
                }
            }
        }
        check_medians();
    }

    NetworKit::index subsets(
            std::vector<h_set> &T, NetworKit::count n, NetworKit::count k, NetworKit::index p) {
        // Stores the subsets of size k of {0, ..., n - 1} in T,
        // starting in position p, and returns p plus their number.
        if (n < k) {
            return p;
        }
        if (k == 0) {
            T[p] = 0;
            return p + 1;
        }
        NetworKit::index q = subsets(T, n - 1, k - 1, p);
        for (auto i = p; i < q; i++) {
            T[i] |= 1 << (n - 1);
        }
        return subsets(T, n - 1, k, q);
    }

    void check_medians() {
        for (std::size_t i = 0; i < median.size(); i++) {
            if (i == 0) {
                assert(median[i] == 0);
                continue;
            }
            std::vector<NetworKit::index> elements;
            auto hset = i;
            for (NetworKit::index j = 0; (1 << j) <= hset; j++) {
                if ((1 << j) & hset) {
                    elements.push_back(j);
                }
            }
            assert(elements.size() > 0);
            assert(elements.at(elements.size() / 2) == median[i]);
        }
    }

    void visit(
            NetworKit::node v, h_set S, const std::vector<NetworKit::node>& upper,
            std::vector<NetworKit::node>& answer) {
        P[depth[v]] = v;  // push current node on stack
        // sup{j' \in down(Dv, Su) : w(Pv(j')) > w(v)}
        auto k = binary_search(P, getWeight(v), down(D[v], S));
        // BUG in the paper: S = down(D[v], S & (1 << (k + 1) - 1) | (1 << depth[v]));
        S = down(D[v], S & ((1 << (k + 1)) - 1) | (1 << depth[v]));
        for (auto i = L[v]; i != NetworKit::none; i = L_next[i]) {
            answer[i] = P[median[down(1 << depth[upper[i]], S)]];
        }
        getTree().forNeighborsOf(v, [&](NetworKit::node child) {
            visit(child, S, upper, answer);
        });
    }

    inline h_set down(const h_set &A, const h_set &B) {
        // Returns A "downarrow" B
        return B & (~(A | B) ^ (A + (A | ~B)));
    }

    // when called, S is "S of `parent(v)`" or \emptyset for root.
    NetworKit::count binary_search(
            const std::vector<NetworKit::node> &P, NetworKit::edgeweight w, int S) {
        // Returns max({j in S | weight[P[j]]>w} union {0})
        // needed for Sv definition on the bottom of paper's 183 page.
        if (S == 0) {
            return 0;
        }
        auto j = median[S];
        // `while |S| > 1` (or, to be more specific, `while S != {j}`).
        for (; S != (1 << j); j = median[S]) {
            S &= (getWeight(P[j]) > w) ? (~((1 << j) - 1)) : ((1 << j) - 1);
        }
        return getWeight(P[j]) > w ? j : 0;
    }
};

namespace Koala {

std::vector<NetworKit::node> get_parent_vector(const NetworKit::Graph &tree) {
    std::vector<NetworKit::node> parent(tree.upperNodeIdBound(), NetworKit::none);
    tree.forEdges([&](NetworKit::node u, NetworKit::node v) {
        parent[v] = u;
    });
    return parent;
}

MinimumSpanningTree::MinimumSpanningTree(
        NetworKit::Graph &graph) : graph(graph), tree(NetworKit::GraphTools::copyNodes(graph)) {}

const NetworKit::Graph& MinimumSpanningTree::getForest() const {
    assureFinished();
    return tree;
}

void MinimumSpanningTree::initialize() {
    tree = NetworKit::GraphTools::copyNodes(graph);
}

void KruskalMinimumSpanningTree::run() {
    GraphTools::assureUndirectedGraph(graph);
    initialize();
    hasRun = true;
    std::vector<NetworKit::WeightedEdge> sorted_edges(
        graph.edgeWeightRange().begin(), graph.edgeWeightRange().end());
    Aux::Parallel::sort(sorted_edges.begin(), sorted_edges.end());
    NetworKit::UnionFind union_find(graph.upperNodeIdBound());
    for (const auto &e : sorted_edges) {
        if (union_find.find(e.u) != union_find.find(e.v)) {
            tree.addEdge(e.u, e.v, e.weight);
            union_find.merge(e.u, e.v);
        }
    }
}

void PrimMinimumSpanningTree::run() {
    GraphTools::assureUndirectedGraph(graph);
    initialize();
    hasRun = true;
    Heap<std::pair<NetworKit::edgeweight, NetworKit::node>> queue;
    queue.push(std::make_pair(0, *(graph.nodeRange().begin())));
    std::unordered_map<NetworKit::node, NetworKit::WeightedEdge> previous;
    while (!queue.empty()) {
        auto v = queue.top().second;
        queue.pop();
        if (!tree.isIsolated(v)) {
            continue;
        }
        const auto &e = previous.find(v);
        if (e != previous.end()) {
            tree.addEdge(e->second.u, e->second.v, e->second.weight);
        }
        graph.forNeighborsOf(v, [&](NetworKit::node u, NetworKit::edgeweight weight) {
            if (tree.isIsolated(u)) {
                queue.push(std::make_pair(-weight, u));
                if (!previous.count(u) || previous[u].weight > weight) {
                    previous[u] = NetworKit::WeightedEdge(u, v, weight);
                }
            }
        });
    }
}

void BoruvkaMinimumSpanningTree::run() {
    GraphTools::assureUndirectedGraph(graph);
    initialize();
    hasRun = true;
    NetworKit::UnionFind union_find(graph.upperNodeIdBound());
    NetworKit::Graph G(graph);
    EdgeMap E;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({NetworKit::Edge(u, v, true), {u, v}});
    });
    iterate(G, tree, union_find, E, std::numeric_limits<NetworKit::count>::max(), false);
}

std::optional<NetworKit::Graph> BoruvkaMinimumSpanningTree::iterate(
        NetworKit::Graph &G, NetworKit::Graph &F,
        NetworKit::UnionFind &union_find, EdgeMap &E,
        NetworKit::count steps, bool get_branching_tree) {
    BranchingTree B;
    if (get_branching_tree) {
        B.initialize(G);
    }
    while (G.numberOfNodes() > 1 && G.numberOfEdges() > 0 && steps-- > 0) {
        G.forNodes([&](NetworKit::node x) {
            const auto &[y, w] = *std::min_element(
                G.weightNeighborRange(x).begin(), G.weightNeighborRange(x).end(),
                [](const auto &e1, const auto &e2) { return e1.second < e2.second; });
            if (get_branching_tree) {
                B.addEdge(x, w);
            }
            NetworKit::node u_prim = union_find.find(x), v_prim = union_find.find(y);
            if (u_prim == v_prim) {
                return;
            }
            union_find.merge(u_prim, v_prim);
            const auto &[u, v] = E[NetworKit::Edge(x, y, true)];
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
                second_pass[union_find.find(e.u)].emplace_back(e);
            }
        }
        auto G_prim = NetworKit::GraphTools::copyNodes(G);
        for (const auto &v : G.nodeRange()) {
            if (union_find.find(v) != v) {
                G_prim.removeNode(v);
            }
        }
        if (get_branching_tree) {
            B.update(G_prim, union_find);
        }
        // It would be nice to clear the E map of the edges that actually don't exist in here
        // of couse E works when you give it proper vertices, but due to not removing edges:
        // - E[u, v] for something that's not a proper edge works and is some garbage
        //      and UB is horrible as there may be an error that uses it...
        // - E is much larger then actually needed... => the complexity may rise for no reason.
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
                assert(u_prim == union_find.find(e.u) && v_prim == union_find.find(e.v));
                E[NetworKit::Edge(u_prim, v_prim, true)] = E[NetworKit::Edge(e.u, e.v, true)];
                left = right;
            }
        }
        std::swap(G, G_prim);
    }
    return get_branching_tree ? std::make_optional(std::move(B.getTree())) : std::nullopt;
}

void KargerKleinTarjanMinimumSpanningTree::run() {
    GraphTools::assureUndirectedGraph(graph);
    initialize();
    hasRun = true;
    NetworKit::Graph G(graph);
    recurse(G, tree);
}

void KargerKleinTarjanMinimumSpanningTree::recurse(NetworKit::Graph &G, NetworKit::Graph &F) {
    NetworKit::UnionFind union_find(G.upperNodeIdBound());
    EdgeMap E;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({NetworKit::Edge(u, v, true), {u, v}});
    });
    while (true) {
        iterate(G, F, union_find, E, 2, false);
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
    for (NetworKit::count i = 1; i < connected_components.numberOfComponents(); i++) {
        subgraph.addEdge(
            components[0][0], components[i][0], std::numeric_limits<NetworKit::edgeweight>::max());
    }
}

void KargerKleinTarjanMinimumSpanningTree::remove_heavy_edges(
        NetworKit::Graph &G, NetworKit::Graph &subforest) {
    NetworKit::UnionFind union_find(G.upperNodeIdBound());
    EdgeMap E;
    subforest.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({NetworKit::Edge(u, v, true), {u, v}});
    });
    NetworKit::Graph branching_tree = *iterate(
        subforest, subforest, union_find, E, std::numeric_limits<NetworKit::count>::max(), true);
    auto branching_tree_augmented = AugmentedGraph(branching_tree);
    Koala::LCA lca(get_parent_vector(branching_tree), branching_tree_augmented.getRoot());
    std::vector<NetworKit::node> upper, lower;
    std::vector<std::tuple<NetworKit::node, NetworKit::node, NetworKit::edgeweight>> edges;
    G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        auto uv = lca.query(u, v);
        upper.push_back(uv), upper.push_back(uv), lower.push_back(u), lower.push_back(v);
        edges.push_back({u, v, w});
    });
    auto answers = branching_tree_augmented.getTreePathMaxima(lower, upper);
    for (NetworKit::index i = 0; i < answers.size(); i += 2) {
        const auto &[u, v, w] = edges[i / 2];
        auto max_first_half = branching_tree_augmented.getWeight(answers[i]);
        auto max_second_half = branching_tree_augmented.getWeight(answers[i + 1]);
        if (w > max_first_half && w > max_second_half) {
            G.removeEdge(u, v);
        }
    }
}

float ChazelleRubinfeldTrevisanMinimumSpanningTree::calculate_approximate_degree(float eps) const {
    const int C = 100;
    NetworKit::count d_est = 0;
    std::uniform_int_distribution<int> uniform_distribution(0, graph.numberOfNodes() - 1);

    for (int i = 0; i < C / eps; ++i) {
        NetworKit::count v = uniform_distribution(generator);
        d_est = std::max(d_est, graph.degreeOut(v));
    }
    return d_est;
}

float ChazelleRubinfeldTrevisanMinimumSpanningTree::calculate_approximate_ccs_count(
        float eps, NetworKit::count bfs_bound, unsigned int w_bound) const {
    NetworKit::count n = graph.numberOfNodes();
    int r = 1/eps/eps + 1;
    float estimate = 0;

    float d_est = calculate_approximate_degree(eps);
    std::uniform_int_distribution<int> uniform_distribution(0, n-1);

    for (int i = 0; i < r; ++i) {
        NetworKit::node u = uniform_distribution(generator);
        float coin_flips = 0;
        float b = 0;
        NetworKit::count max_d = 0;
        NetworKit::count visited_edges = 0;
        std::unordered_set<NetworKit::node> vis;
        std::queue<NetworKit::node> q;

        if (graph.degreeOut(u) > d_est) {
            continue;
        }

        NetworKit::count d_actual = 0;

        graph.forEdgesOf(u, [w_bound, &d_actual](
                NetworKit::node, NetworKit::node, NetworKit::edgeweight ew, NetworKit::edgeid) {
            if (ew > w_bound) {
                return;
            }
            d_actual += 1;
        });

        // isolated vertex
        if (d_actual == 0) {
            // self loop
            estimate += 2;
            continue;
        }


        vis.insert(u);
        graph.forEdgesOf(u, [w_bound, &visited_edges, &q](
                NetworKit::node, NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid) {
            if (ew > w_bound) {
                return;
            }
            q.push(v);
            visited_edges += 1;
        });

        while (vis.size() < bfs_bound && max_d <= d_est) {
            coin_flips += 1;
            if (distribution(generator)) {
                break;
            }

            auto prev_visited_edges = 2 * visited_edges;
            while (!q.empty() && visited_edges <= prev_visited_edges) {
                auto v = q.front();
                q.pop();

                if (vis.contains(v)) {
                    continue;
                }
                vis.insert(v);

                graph.forEdgesOf(v, [w_bound, &visited_edges, &q, &vis](
                        NetworKit::node, NetworKit::node v, NetworKit::edgeweight ew,
                        NetworKit::edgeid) {
                    if (ew > w_bound) {
                        return;
                    }
                    if (vis.contains(v)) {
                        return;
                    }
                    q.push(v);
                    visited_edges += 1;
                });
            }
            if (q.empty()) {
                if (visited_edges == 0) {
                    b = 2;
                } else {
                    b = powl(2, coin_flips) * d_actual / visited_edges;
                }
                break;
            }
        }

        estimate += b;
    }

    return estimate * n / 2.0f / r;
}

NetworKit::edgeweight
ChazelleRubinfeldTrevisanMinimumSpanningTree::calculate_approximate_tree_weight(
        float eps, unsigned int w) const {
    NetworKit::edgeweight tree_weight = static_cast<NetworKit::edgeweight>(
        graph.numberOfNodes()) - static_cast<NetworKit::edgeweight>(w);

    for (unsigned int w_bound = 1; w_bound < w; ++w_bound) {
        tree_weight += calculate_approximate_ccs_count(eps, 4 / eps, w_bound);
    }
    return tree_weight;
}

void ChazelleRubinfeldTrevisanMinimumSpanningTree::run() {
    throw std::runtime_error("this function is a stub, use run with parameters instead.");
}

void ChazelleRubinfeldTrevisanMinimumSpanningTree::run(unsigned int w, float eps) {
    GraphTools::assureUndirectedGraph(graph);
    assert(0 < eps && eps < 0.5);
    assert(0 < w);

    tree_weight = calculate_approximate_tree_weight(eps, w);
}

const NetworKit::Graph& ChazelleRubinfeldTrevisanMinimumSpanningTree::getForest() const {
    throw std::runtime_error(
        "getForest method not supported for approximate MST algorithm. Try getTreeWeight");
}

NetworKit::edgeweight ChazelleRubinfeldTrevisanMinimumSpanningTree::getTreeWeight() const {
    return tree_weight;
}

void MinimumSpanningTree::check() const {
    assureFinished();
    assert(tree.numberOfNodes() == tree.numberOfEdges() + 1);
    auto connected_components = NetworKit::ConnectedComponents(tree);
    connected_components.run();
    assert(connected_components.getComponents().size() == 1);

    NetworKit::Graph tree_copy(tree);
    NetworKit::UnionFind union_find(graph.upperNodeIdBound());
    EdgeMap E;
    tree_copy.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({NetworKit::Edge(u, v, true), {u, v}});
    });
    // Note: Boruvka runs in linear time on trees
    NetworKit::Graph branching_tree = *BoruvkaMinimumSpanningTree::iterate(
        tree_copy, tree_copy, union_find, E, std::numeric_limits<NetworKit::count>::max(), true);
    std::vector<std::tuple<NetworKit::node, NetworKit::node, NetworKit::edgeweight>> G_minus_M;
    G_minus_M.reserve(graph.numberOfEdges() - tree.numberOfEdges());
    graph.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (!tree.hasEdge(u, v)) {
            G_minus_M.push_back({u, v, w});
        }
    });
    auto branching_tree_augmented = AugmentedGraph(branching_tree);
    Koala::LCA lca(get_parent_vector(branching_tree), branching_tree_augmented.getRoot());
    std::vector<NetworKit::node> lower, upper;
    lower.reserve(2 * G_minus_M.size()), upper.reserve(2 * G_minus_M.size());
    for (const auto &[u, v, w] : G_minus_M) {
        auto uv = lca.query(u, v);
        lower.push_back(u), lower.push_back(v), upper.push_back(uv), upper.push_back(uv);
    }
    auto answers = branching_tree_augmented.getTreePathMaxima(lower, upper);
    for (size_t i = 0; i < answers.size(); i++) {
        assert(branching_tree_augmented.getWeight(answers[i]) <= std::get<2>(G_minus_M[i / 2]));
    }
}

void Chazelle2000MinimumSpanningTree::run() {
    GraphTools::assureUndirectedGraph(graph);
    initialize();
    auto G = graph;
    tree = mst(G, 10);
    hasRun = true;
}


/**
 * A wrapper for BoruvkaMinimumSpanningTree::iterate function
 * that's used in the Chazelle (2000) algorithm implementation.
 *
 * The function applies c consecutive Boruvka steps
 * returns a tuple [G_minor, edge_id, forest]
 * G_minor - a minor of the input graph G - the result of edge contraction.
 *
 * G_minor_to_G - a map that maps G_minor edges to edges in G
 *          [u, v] -> [uG, vG]
 *
 * forest - edges from the input graph G which were contracted.
 */
std::tuple<NetworKit::Graph, EdgeMap, NetworKit::Graph>
    Chazelle2000MinimumSpanningTree::boruvka_steps(NetworKit::Graph G, int c) {
    assert(c >= 1);

    auto F = NetworKit::GraphTools::copyNodes(G);
    auto uf = NetworKit::UnionFind(graph.upperNodeIdBound());
    EdgeMap G_minor_to_G;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        G_minor_to_G.insert({NetworKit::Edge(u, v, true), {u, v}});
    });
    auto G_minor = G;
    BoruvkaMinimumSpanningTree::iterate(G_minor, F, uf, G_minor_to_G, c, false);

    auto map1 = NetworKit::GraphTools::getContinuousNodeIds(G_minor);
    auto G_minor_compact = NetworKit::GraphTools::getCompactedGraph(G_minor, map1);
    EdgeMap G_minor_remapped_to_G;

    for (auto const &[key, val] : G_minor_to_G) {
        auto old_u = key.u;
        auto old_v = key.v;
        auto new_u_it = map1.find(old_u);
        auto new_v_it = map1.find(old_v);

        if (new_u_it == map1.end() || new_v_it == map1.end()) {
            continue;
        }

        auto new_u = new_u_it->second;
        auto new_v = new_v_it->second;

        G_minor_remapped_to_G[NetworKit::Edge(new_u, new_v, true)] = val;
    }
    return {G_minor_compact, G_minor_remapped_to_G, F};
}

namespace {
struct edge {
    int u;          // G0 node
    int v;          // G0 node
    NetworKit::edgeweight key;        //
    NetworKit::edgeweight ckey;       // for soft heap
    int id;         // for soft heap delete
    bool removed;   // lazy removal from SoftHeap
    bool corrupted;

    friend bool operator<(edge e0, edge e1) {
        if (e0.key == e1.key) return e0.id < e1.id;
        return e0.key < e1.key;
    }

    std::string to_string() const {
        std::ostringstream oss;
        oss << "edge{"
            << "u=" << u
            << ", v=" << v
            << ", key=" << key
            << ", ckey=" << ckey
            << ", id=" << id
            << ", removed=" << std::boolalpha << removed
            << ", corrupted=" << std::boolalpha << corrupted
            << "}";
        return oss.str();
    }
};
}  // namespace

// Max number of nodes ~ billion
std::pair<int, std::vector<int>> t_hierarchy_size(int n, int) {
    if (n < 3) {
        return {{3}, {10, 4, 1}};
    }
    // Arbitrary sequence for now that grows really fast -- no point in using Ackermann sequence.
    int desired_number_of_leaves[7] = {
        1, 1 << 2, 1 << 6, 1 << 16, 1 << 30, 1 << 30, 1 << 30};
    int desired_number_of_children[7];
    desired_number_of_children[0] = 1;
    for (int i = 1; i < 7; ++i) {
        desired_number_of_children[i] =
            desired_number_of_leaves[i] / desired_number_of_leaves[i - 1];
    }
    std::vector<int> desired_size;
    int d = 0;
    for (; desired_number_of_leaves[d] < n; ++d) {}
    for (int i = d; i >= 0; --i) desired_size.push_back(desired_number_of_children[i]);

    return {d, desired_size};
}

NetworKit::Graph Chazelle2000MinimumSpanningTree::msf(NetworKit::Graph G, int t) {
    auto tree = NetworKit::GraphTools::copyNodes(G);
    auto connected_components = NetworKit::ConnectedComponents(G);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (NetworKit::count i = 0; i < connected_components.numberOfComponents(); i++) {
        auto subgraph = NetworKit::GraphTools::subgraphFromNodes(G, [](const auto &components) {
            std::unordered_set<NetworKit::node> s;
            for (auto n : components) {
                s.insert(n);
            }
            return s;
        }(components[i]));

        std::unordered_map<NetworKit::node, NetworKit::node> mapping_to;
        std::unordered_map<NetworKit::node, NetworKit::node> mapping_from;
        subgraph.forNodes([&](NetworKit::node u) {
            NetworKit::node n = mapping_to.size();
            mapping_to[u] = n;
            mapping_from[n] = u;
        });
        subgraph = NetworKit::GraphTools::getCompactedGraph(subgraph, mapping_to);
        auto tree1 = mst(subgraph, t);
        tree1.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew) {
            tree.addEdge(mapping_from[u], mapping_from[v], ew);
        });
    }
    return tree;
}

NetworKit::Graph Chazelle2000MinimumSpanningTree::mst(NetworKit::Graph G, int t) {
    // The paper bounds the number of bad edges by 8m'/c + d^3n' and requires
    // it to be at most m'/2 + d^3n'. Its smallest hierarchy target is S(t, 1)^3 = 8.
    static constexpr int C = 16;
    static constexpr int MIN_NUMBER_NODES = 8;

    // Expects that G is connected and undirected
    assert(Koala::GraphTools::isConnected(G));
    assert(!G.isDirected());
    assert(!Koala::GraphTools::hasMultiEdges(G));

    // [STEP 1]
    if (t <= 1 || G.numberOfNodes() < MIN_NUMBER_NODES) {
        BoruvkaMinimumSpanningTree boruvka(G);
        boruvka.run();
        return boruvka.getForest();
    }

    // [STEP 2]
    auto [G0, edge_G0_to_G, forest] = boruvka_steps(G, C);
    if (G0.numberOfEdges() == 0) {
        return forest;
    }

    edge dummy_edge{
        0, 0, std::numeric_limits<NetworKit::edgeweight>::max(),
        std::numeric_limits<NetworKit::edgeweight>::max(), 0, true, true};
    std::vector<std::vector<SoftHeap<edge*>>> heaps;
    std::unordered_set<NetworKit::Edge> contracted_edges;
    std::unordered_set<NetworKit::Edge> bad_edges;
    std::vector<edge> edges(G0.numberOfEdges());
    std::vector<std::set<int>> Cz;
    std::vector<int> parent(G0.numberOfNodes());
    std::vector<bool> fusion_node(G0.numberOfNodes(), false);
    std::vector<std::vector<std::pair<NetworKit::Edge, NetworKit::edgeweight>>> min_link;
    std::unordered_map<NetworKit::Edge, int> G0_edge_id;
    std::vector<bool> visited(G0.numberOfNodes(), false);
    std::vector<std::stack<int>> border(G0.numberOfNodes(), std::stack<int>());
    assert(!G.isDirected());
    assert(!G0.isDirected());
    auto [d, desired_size] = t_hierarchy_size(G0.numberOfNodes(), G0.numberOfEdges());


    // we have Cz0, Cz1, .... Czk
    auto k = [&]() { return Cz.size() - 1; };

    auto min_border_edge = [&]() {
        NetworKit::edgeweight min_key = std::numeric_limits<NetworKit::edgeweight>::max();
        int mini = -1, minj = -1;

        for (std::size_t i = 0; i < heaps.size(); ++i) {
            for (std::size_t j = 0; j < heaps[i].size(); ++j) {
                auto ret = heaps[i][j].lookupMin();
                if (std::holds_alternative<bool>(ret)) continue;
                auto val = std::get<edge*>(ret);
                if (val->ckey < min_key) {
                    mini = i;
                    minj = j;
                    min_key = val->ckey;
                }
            }
        }
        edge* min_edge = heaps[mini][minj].extractMin();
        return std::tuple{*min_edge, 0, 0};
    };

    auto find_and_delete_edge_from_heaps = [&](NetworKit::node u, NetworKit::node v) {
        auto eid = G0_edge_id[NetworKit::Edge(u, v, true)];
        edges[eid].removed = true;
        if (edges[eid].corrupted) {
            bad_edges.insert(NetworKit::Edge(u, v, true));
        }
        return edges[eid];
    };

    auto node_top_parent = [&](NetworKit::node u) {
        int p = static_cast<int>(u);
        while (parent[p] != p) p = parent[p];
        return p;
    };

    auto czi_of_node = [&](NetworKit::node u) {
        int p = node_top_parent(u);
        for (std::size_t i = 0; i < Cz.size(); ++i) {
            if (Cz[i].contains(p)) {
                return i;
            }
        }
        throw std::runtime_error("u node should be found");
    };

    auto insert_new_border_edge = [&](edge* e, NetworKit::node u, NetworKit::node v) {
        int j = Cz.size();
        if (border[v].empty()) {
            heaps[0][j] = std::move(insert(std::move(heaps[0][j]), e));
        } else {
            int i = 1 + czi_of_node(border[v].top());
            heaps[i][j] = std::move(insert(std::move(heaps[i][j]), e));
        }
        border[v].push(u);
    };

    auto leftmost_smaller_min_link = [&](int ckey) {
        NetworKit::Edge ab;

        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = i + 1; j <= k(); ++j) {
                if (min_link[i][j].second < 0) continue;
                if (min_link[i][j].second <= ckey) {
                    return std::tuple<int, int, NetworKit::Edge>{i, j, min_link[i][j].first};
                }
            }
        }
        return std::tuple<int, int, NetworKit::Edge>{-1, -1, ab};
    };

    std::set<edge> inside_edges;

    auto update_min_links = [&]() {
        int last_i = k();
        for (int i = 0; i < last_i; ++i) {
            int best_j = last_i;
            NetworKit::edgeweight best_val = min_link[i][last_i].second;
            for (size_t j = last_i + 1; j < min_link[i].size(); ++j) {
                if (min_link[i][j].second < 0) continue;
                if (min_link[i][j].second >= best_val && best_val >= 0) continue;
                best_j = j;
                best_val = min_link[i][j].second;
            }
            min_link[i][last_i] = min_link[i][best_j];
        }

        for (std::size_t i = 0; i < d; ++i) {
            for (std::size_t j = 0; j < d; ++j) {
                if (i >= k() || j > k()) {
                    min_link[i][j] = {
                        NetworKit::Edge(), std::numeric_limits<NetworKit::edgeweight>::max()};
                }
            }
        }
    };

    auto pop_heaps = [&]() {
        int k = Cz.size();
        heaps[k - 1][k - 1] = std::move(meld(std::move(heaps[k - 1][k - 1]),
            std::move(heaps[k - 1][k])));
        heaps[k - 1][k - 1] = std::move(meld(std::move(heaps[k - 1][k - 1]),
            std::move(heaps[k][k])));

        for (int i = 0; i < k - 1; ++i) {
            heaps[i][k - 1] = std::move(meld(std::move(heaps[i][k - 1]), std::move(heaps[i][k])));
        }
        heaps.pop_back();
        for (auto &h : heaps) {
            h.pop_back();
        }
    };

    auto retraction = [&]() {
        int new_p = parent.size();
        parent.push_back(new_p);
        fusion_node.push_back(false);
        auto vertices = Cz[Cz.size() - 1];
        pop_heaps();
        Cz.pop_back();
        if (vertices.size() == 1) {
            parent.pop_back();
            fusion_node.pop_back();
            new_p = *vertices.begin();
        } else {
            for (auto v : vertices) {
                parent[v] = new_p;
            }
        }
        Cz[Cz.size() - 1].insert(new_p);
        update_min_links();
    };

    auto fusion = [&]() -> edge {
        auto [uv_edge, mini, minj] = min_border_edge();
        NetworKit::node uG0{uv_edge.u}, vG0{uv_edge.v};
        auto [link_i, link_j, ab] = leftmost_smaller_min_link(uv_edge.ckey);
        if (link_i == -1) return uv_edge;
        auto [a, b] = ab;
        auto ap = node_top_parent(a);
        auto bp = node_top_parent(b);
        if (!Cz[link_i].contains(ap)) {
            std::swap(a, b);
            std::swap(ap, bp);
        }

        size_t new_p = parent.size();
        parent.push_back(new_p);
        fusion_node.push_back(true);
        while (k() > link_i) {
            for (auto v : Cz[k()]) {
                parent[v] = new_p;
            }
            pop_heaps();
            Cz.pop_back();
        }
        Cz[k()].erase(ap);
        parent[ap] = new_p;
        fusion_node[new_p] = true;
        Cz[k()].insert(new_p);

        update_min_links();
        return uv_edge;
    };

    auto update_min_links_with_edge = [&](edge e) {
        int i = czi_of_node(e.u), j = czi_of_node(e.v);
        if (i > j) std::swap(i, j);
        if (e.ckey < min_link[i][j].second) {
            min_link[i][j] = {NetworKit::Edge(e.u, e.v, true), e.ckey};
        }
    };

    int extensions = 0;
    auto extension = [&](edge e) {
        extensions += 1;
        NetworKit::node u = e.u, v = e.v;
        if (!visited[u]) std::swap(u, v);

        heaps.push_back(std::vector<SoftHeap<edge*>>());
        for (std::size_t i = 0; i <= Cz.size(); ++i) heaps[heaps.size() - 1].push_back(
            SoftHeap<edge*>(&dummy_edge, 0.1));
        Cz.push_back({v});
        for (std::size_t i = 0; i < heaps.size(); ++i) {
            heaps[i].push_back(SoftHeap<edge*>(&dummy_edge, 0.1));
        }
        visited[v] = true;
        G0.forEdgesOf(v, [&](NetworKit::node, NetworKit::node w) {
            if (visited[w]) {
                edge e = find_and_delete_edge_from_heaps(v, w);
                update_min_links_with_edge(e);
                inside_edges.insert(e);
            } else {
                int eid = G0_edge_id[NetworKit::Edge(v, w, true)];
                insert_new_border_edge(&edges[eid], v, w);
            }
        });
        update_min_links();
    };

    auto should_retract = [&]() {
        return Cz.size() == d || (Cz[k()].size() >= desired_size[k()] && k() > 0);
    };

    auto should_finish = [&]() {
        return extensions >= G0.numberOfNodes() - 1;
    };

    auto initialization = [&]() {
        int eid_now = 0;
        G0.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew) {
            edges[eid_now] = edge{u, v, ew, ew, eid_now, false, false};
            G0_edge_id[NetworKit::Edge(u, v, true)] = eid_now;
            eid_now += 1;
        });

        for (size_t i = 0; i < parent.size(); i++) parent[i] = i;
        Cz.push_back({0});
        heaps.push_back({SoftHeap<edge*>(&dummy_edge, 0.1), SoftHeap<edge*>(&dummy_edge, 0.1)});
        heaps.push_back({SoftHeap<edge*>(&dummy_edge, 0.1), SoftHeap<edge*>(&dummy_edge, 0.1)});
        G0.forEdgesOf(0, [&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew) {
            int eid = G0_edge_id[NetworKit::Edge(u, v, true)];
            insert_new_border_edge(&edges[eid], u, v);
        });
        min_link = std::vector<
            std::vector<std::pair<NetworKit::Edge, NetworKit::edgeweight>>>(
                d, std::vector<std::pair<NetworKit::Edge, NetworKit::edgeweight>>(
                       d, {NetworKit::Edge(),
                           std::numeric_limits<NetworKit::edgeweight>::max()}));
        visited[0] = true;
    };

    initialization();
    while (!should_finish()) {
        if (should_retract()) {
            retraction();
        } else {
            edge e = fusion();
            extension(e);
        }
    }

    while (Cz.size() > 1) {
        size_t new_p = parent.size();
        parent.push_back(new_p);
        fusion_node.push_back(false);
        auto vertices = Cz[Cz.size() - 1];

        heaps[k() - 1][k() - 1] = std::move(
            meld(std::move(heaps[k()][k()]), std::move(heaps[k() - 1][k() - 1])));
        heaps.pop_back();
        Cz.pop_back();

        for (auto v : vertices) {
            parent[v] = new_p;
        }
        Cz[Cz.size() - 1].insert(new_p);
    }
    size_t root = parent.size();
    parent.push_back(root);
    fusion_node.push_back(false);
    for (int v : Cz[0]) {
        parent[v] = root;
    }

    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (!edges[i].corrupted) continue;
        bad_edges.insert(NetworKit::Edge(edges[i].u, edges[i].v, true));
    }

    // [STEP 4]
    std::vector<std::map<NetworKit::node, NetworKit::node>> from_G0_maps(parent.size());
    std::vector<std::unordered_map<NetworKit::Edge, int>> to_G0_edge_id(parent.size());
    std::vector<int> depth(parent.size(), 0);

    for (NetworKit::count i = 0; i < G0.numberOfNodes(); ++i) {
        int v = i;
        while (v != parent[v]) {
            v = parent[v];
            depth[i] += 1;
        }

        int d = depth[i];
        v = i;
        while (v != parent[v]) {
            int prev_v = v;
            v = parent[v];
            d -= 1;

            depth[v] = d;
            int new_v = from_G0_maps[v].size();
            from_G0_maps[v].insert({prev_v, new_v});
        }
    }

    std::vector<NetworKit::Graph> Cz_graphs(parent.size(), NetworKit::Graph(0, true));
    for (NetworKit::count i = G0.numberOfNodes(); i < parent.size(); ++i) {
        Cz_graphs[i] = NetworKit::Graph(from_G0_maps[i].size(), true);
    }

    auto lca_lr = [&](int l, int r) {
        int prev_l = l, prev_r = r;
        while (depth[l] > depth[r]) {
            prev_l = l;
            l = parent[l];
        }
        while (depth[l] < depth[r]) {
            prev_r = r;
            r = parent[r];
        }

        while (l != r) {
            prev_l = l; prev_r = r;
            l = parent[l]; r = parent[r];
        }

        return std::tuple<int, int, int>{l, prev_l, prev_r};
    };

    std::vector<std::unordered_map<NetworKit::Edge, NetworKit::edgeweight>> smallest_edge(
        Cz_graphs.size());

    G0.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew) {
        if (bad_edges.contains(NetworKit::Edge(u, v, true))) return;
        auto [lca, l, r] = lca_lr(u, v);
        auto minor_l = from_G0_maps[lca][l];
        auto minor_r = from_G0_maps[lca][r];
        auto np = NetworKit::Edge(minor_l, minor_r, true);
        if (!smallest_edge[lca].contains(np)) {
            smallest_edge[lca][np] = ew;
            to_G0_edge_id[lca][np] = G0_edge_id[NetworKit::Edge(u, v, true)];
        }
        if (smallest_edge[lca][np] > ew) {
            smallest_edge[lca][np] = ew;
            to_G0_edge_id[lca][np] = G0_edge_id[NetworKit::Edge(u, v, true)];
        }
    });
    for (NetworKit::count i = G0.numberOfNodes(); i < Cz_graphs.size(); ++i) {
        for (auto [key, val] : smallest_edge[i]) {
            Cz_graphs[i].addEdge(key.u, key.v, val);
        }
    }

    NetworKit::Graph G0_copy = G0;
    BoruvkaMinimumSpanningTree bbb(G0_copy);
    bbb.run();
    auto actual_mst = bbb.getForest();
    std::unordered_set<NetworKit::Edge> actual_mst_edges;
    actual_mst.forEdges([&](NetworKit::node u, NetworKit::node v) {
        actual_mst_edges.insert(NetworKit::Edge(u, v, true));
    });

    for (NetworKit::count i = G0.numberOfNodes(); i < Cz_graphs.size(); ++i) {
        if (fusion_node[i]) continue;
        Cz_graphs[i] = msf(Cz_graphs[i], t - 1);
    }

    auto F = NetworKit::GraphTools::copyNodes(G0);
    for (NetworKit::count i = G0.numberOfNodes(); i < Cz_graphs.size(); ++i) {
        Cz_graphs[i].forEdges([&](NetworKit::node u, NetworKit::node v) {
            const auto &e = edges[to_G0_edge_id[i][NetworKit::Edge(u, v, true)]];
            F.addEdge(e.u, e.v, e.key);
        });
    }

    for (auto [x, y] : bad_edges) {
        F.addEdge(x, y, edges[G0_edge_id[NetworKit::Edge(x, y, true)]].key);
    }

    // [STEP 5]
    NetworKit::Graph res = mst(F, t);

    res.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew) {
        auto [uG, vG] = edge_G0_to_G[NetworKit::Edge(u, v, true)];
        forest.addEdge(uG, vG, ew);
    });

    for (size_t i = 0; i < parent.size(); ++i) {
        int v = i;
        while (v != parent[v]) {
            v = parent[v];
        }
    }

    return forest;
}

}  // namespace Koala
