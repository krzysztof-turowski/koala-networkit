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
#include <unordered_set>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <limits>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <stack>
#include <variant>
#include <string>
#include <tuple>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <structures/Heap.hpp>
#include <structures/LCA.hpp>
// #include <structures/heap/SoftHeap.hpp>

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
    // `L` and `Lnext` combined form a linkedlist of "indices of queries" (query_ids).
    std::vector<NetworKit::index> L, Lnext;

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
        L.resize(n, NetworKit::none), Lnext.resize(upper.size(), NetworKit::none);
        for (auto i = 0; i < lower.size(); i++) {  // distribute queries to lower nodes.
            // L[u] - beginning of a linked list query_ids that have lower node `u`.
            // L[u] == index of the first query. Lnext[L[u]] - index of the next query.
            // L[] is indexed by nodes and Lnext[] is indexed by query_ids.
            // L[u] points to the beginning of the list for `u`. This list is stored in `Lnext[]`.
            Lnext[i] = L[lower[i]], L[lower[i]] = i;
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
            const std::vector<NetworKit::node>& upper) {
        depth[u] = u_depth;
        if (u_depth > height) {
            height = u_depth;
        }
        for (auto i = L[u]; i != NetworKit::none; i = Lnext[i]) {
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
        for (int i = 0; i < median.size(); i++) {
            if (i == 0) {
                assert(median[i] == 0);
                continue;
            }
            std::vector<int> elements;
            auto hset = i;
            for (int j = 0; (1 << j) <= hset; j++) {
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
        int k = binary_search(P, getWeight(v), down(D[v], S));
        // BUG in the paper: S = down(D[v], S & (1 << (k + 1) - 1) | (1 << depth[v]));
        S = down(D[v], S & ((1 << (k + 1)) - 1) | (1 << depth[v]));
        for (auto i = L[v]; i != NetworKit::none; i = Lnext[i]) {
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
    NetworKit::count binary_search(const std::vector<NetworKit::node> &P, double w, int S) {
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
    std::map<NodePair, NodePair> E;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({std::minmax(u, v), {u, v}});
    });
    iterate(G, *tree, union_find, E, std::numeric_limits<NetworKit::count>::max(), false);
}

std::optional<NetworKit::Graph> BoruvkaMinimumSpanningTree::iterate(
        NetworKit::Graph &G, NetworKit::Graph &F,
        NetworKit::UnionFind &union_find, std::map<NodePair, NodePair> &E,
        NetworKit::count steps, bool get_branching_tree) {
    BranchingTree B;
    std::unordered_map<NetworKit::node, NetworKit::node> V_B;
    if (get_branching_tree) {
        B.initialize(G);
    }
    while (G.numberOfNodes() > 1 && G.numberOfEdges() > 0 && steps-- > 0) {
        std::unordered_map<NetworKit::node, NetworKit::edgeweight> B_edges;
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
                E[{u_prim, v_prim}] = E[std::minmax(e.u, e.v)];
                left = right;
            }
        }
        std::swap(G, G_prim);
    }
    return get_branching_tree ? std::make_optional(std::move(B.getTree())) : std::nullopt;
}

void KargerKleinTarjanMinimumSpanningTree::run() {
    hasRun = true;
    NetworKit::Graph G(*graph);
    recurse(G, *tree);
}

void KargerKleinTarjanMinimumSpanningTree::recurse(NetworKit::Graph &G, NetworKit::Graph &F) {
    NetworKit::UnionFind union_find(G.upperNodeIdBound());
    std::map<NodePair, NodePair> E;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({std::minmax(u, v), {u, v}});
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
    for (auto i = 1; i < connected_components.numberOfComponents(); i++) {
        subgraph.addEdge(
            components[0][0], components[i][0], std::numeric_limits<NetworKit::edgeweight>::max());
    }
}

void KargerKleinTarjanMinimumSpanningTree::remove_heavy_edges(
        NetworKit::Graph &G, NetworKit::Graph &subforest) {
    NetworKit::UnionFind union_find(G.upperNodeIdBound());
    std::map<NodePair, NodePair> E;
    subforest.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({std::minmax(u, v), {u, v}});
    });
    NetworKit::Graph branching_tree = *iterate(
        subforest, subforest, union_find, E, std::numeric_limits<NetworKit::count>::max(), true);
    auto branching_tree_augmented = AugmentedGraph(branching_tree);
    Koala::LCA<AugmentedGraph> lca(branching_tree_augmented);
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

float ChazelleRubinfeldTrevisanMinimumSpanningTree::calculateApproximateDegree(float eps) const {
    const int C = 100;
    int d_est = 0;
    std::uniform_int_distribution<int> uniform_distribution(0, graph->numberOfNodes() - 1);

    for (int i = 0; i < C / eps; ++i) {
        int v = uniform_distribution(generator);
        d_est = std::max(d_est, static_cast<int>(graph->degreeOut(v)));
    }
    return d_est;
}

float ChazelleRubinfeldTrevisanMinimumSpanningTree::calculateApproximateCCsCount(
    float eps, int bfs_bound, unsigned int w, unsigned int w_bound
) const {
    int n = graph->numberOfNodes();
    int r = 1/eps/eps + 1;
    float estimate = 0;

    float d_est = calculateApproximateDegree(eps);
    std::uniform_int_distribution<int> uniform_distribution(0, n-1);

    for (int i = 0; i < r; ++i) {
        int u = uniform_distribution(generator);
        float coin_flips = 0;
        float b = 0;
        int max_d = 0;
        int visited_edges = 0;
        std::unordered_set<int> vis;
        std::queue<int> q;

        if (graph->degreeOut(u) > d_est) {
            continue;
        }

        int d_actual = 0;

        graph->forEdgesOf(u, [this, w_bound, &d_actual]
            (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid id) {
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
        graph->forEdgesOf(u, [this, w_bound, &visited_edges, &q](NetworKit::node u,
            NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid id){
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

            int prev_visited_edges = 2*visited_edges;
            while (!q.empty() && visited_edges <= prev_visited_edges) {
                int v = q.front();
                q.pop();

                if (vis.contains(v)) {
                    continue;
                }
                vis.insert(v);

                graph->forEdgesOf(v, [this, w_bound, &visited_edges, &q, &vis](NetworKit::node u,
                    NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid id) {
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
                    b = powl(2, coin_flips) *
                        static_cast<float>(d_actual) / static_cast<float>(visited_edges);
                }
                break;
            }
        }

        estimate += b;
    }

    return estimate * static_cast<float>(n) / 2.0f / r;
}

float ChazelleRubinfeldTrevisanMinimumSpanningTree::calculateApproximateTreeWeight(
    float eps, unsigned int w
) const {
    float approx = graph->numberOfNodes() - static_cast<float>(w);

    for (int w_bound = 1; w_bound < w; ++w_bound) {
        float ccs = calculateApproximateCCsCount(eps, 4/eps, w, w_bound);
        approx += ccs;
    }
    return approx;
}

void ChazelleRubinfeldTrevisanMinimumSpanningTree::run() {
    throw std::runtime_error("this function is a stub, use run with parameters instead.");
}

void ChazelleRubinfeldTrevisanMinimumSpanningTree::run(unsigned int w, float eps) {
    assert(0 < eps && eps < 0.5);
    assert(0 < w);

    treeWeight = calculateApproximateTreeWeight(eps, w);
}

const NetworKit::Graph& ChazelleRubinfeldTrevisanMinimumSpanningTree::getForest() const {
    throw std::runtime_error(
        "getForest method not supported for approximate MST algorithm. Try getTreeWeight");
}

float ChazelleRubinfeldTrevisanMinimumSpanningTree::getTreeWeight() const {
    return treeWeight;
}

void MinimumSpanningTree::check() const {
    assureFinished();
    assert(tree->numberOfNodes() == tree->numberOfEdges() + 1);
    auto connected_components = NetworKit::ConnectedComponents(*tree);
    connected_components.run();
    assert(connected_components.getComponents().size() == 1);

    NetworKit::Graph tree_copy(*tree);
    NetworKit::UnionFind union_find(graph->upperNodeIdBound());
    std::map<NodePair, NodePair> E;
    tree_copy.forEdges([&](NetworKit::node u, NetworKit::node v) {
        E.insert({std::minmax(u, v), {u, v}});
    });
    // Note: Boruvka runs in linear time on trees
    NetworKit::Graph branching_tree = *BoruvkaMinimumSpanningTree::iterate(
        tree_copy, tree_copy, union_find, E, std::numeric_limits<NetworKit::count>::max(), true);
    std::vector<std::tuple<NetworKit::node, NetworKit::node, NetworKit::edgeweight>> G_minus_M;
    G_minus_M.reserve(graph->numberOfEdges() - tree->numberOfEdges());
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
         if (!tree->hasEdge(u, v)) {
             G_minus_M.push_back({u, v, w});
         }
    });
    auto branching_tree_augmented = AugmentedGraph(branching_tree);
    auto lca = Koala::LCA(branching_tree_augmented);
    std::vector<NetworKit::node> lower, upper;
    lower.reserve(2 * G_minus_M.size()), upper.reserve(2 * G_minus_M.size());
    for (const auto &[u, v, w] : G_minus_M) {
        auto uv = lca.query(u, v);
        lower.push_back(u), lower.push_back(v), upper.push_back(uv), upper.push_back(uv);
    }
    auto answers = branching_tree_augmented.getTreePathMaxima(lower, upper);
    for (NetworKit::index i = 0; i < answers.size(); i++) {
        assert(branching_tree_augmented.getWeight(answers[i]) <= std::get<2>(G_minus_M[i / 2]));
    }
}
namespace {
/**
 * Assumes that both graphs have the same nodes
 */
NetworKit::Graph graphSum(const NetworKit::Graph& g1, const NetworKit::Graph& g2) {
    NetworKit::Graph g;
    g1.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
        g.addEdge(u, v, ew);
    });
    g2.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
        g.addEdge(u, v, ew);
    });
    return g;
}
}  // namespace

void Chazelle2000MinimumSpanningTree::run() {
    auto G = graph.value();
    tree = mst(G, 10);
    hasRun = true;
}


/**
 * A wrapper for BoruvkaMinimumSpanningTree::iterate function
 * that's used in the Chazelle(2000) algorithm implementation.
 *
 * The function applies c consecutive Boruvka steps
 * returns a tuple [gMinor, edgeId, forest]
 * gMinor - a minor of the input graph G - the result of edge contraction.
 *
 * gMinor_to_G - a map that maps gMinor edges to edges in G
 *          [u, v] -> [uG, vG]
 * 
 * forest - edges from the input graph G which were contracted.
 */
std::tuple<NetworKit::Graph, std::map<Koala::MinimumSpanningTree::NodePair,
    Koala::MinimumSpanningTree::NodePair>, NetworKit::Graph>
    Chazelle2000MinimumSpanningTree::boruvkaSteps(NetworKit::Graph G, int c) {
    assert(c >= 1);

    auto F = NetworKit::GraphTools::copyNodes(G);
    auto uf = NetworKit::UnionFind(graph->upperNodeIdBound());
    std::map<NodePair, NodePair> gMinor_to_G;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        gMinor_to_G.insert({std::minmax(u, v), {u, v}});
    });
    auto gMinor = G;
    BoruvkaMinimumSpanningTree::iterate(gMinor, F, uf, gMinor_to_G, c, false);
    using NetworKit::node; using NetworKit::edgeweight;

    auto map1 = NetworKit::GraphTools::getContinuousNodeIds(gMinor);
    auto gMinorCompact = NetworKit::GraphTools::getCompactedGraph(gMinor, map1);
    std::map<NodePair, NodePair> gMinorRemapped_to_G;

    for (auto const& [key, val] : gMinor_to_G) {
        auto oldU = key.first;
        auto oldV = key.second;
        auto newU_it = map1.find(oldU);
        auto newV_it = map1.find(oldV);

        if (newU_it == map1.end() || newV_it == map1.end()) {
            continue;
        }

        auto newU = newU_it->second;
        auto newV = newV_it->second;

        gMinorRemapped_to_G[std::minmax(newU, newV)] = val;
    }
    return { gMinorCompact, gMinorRemapped_to_G, F };
}

namespace {
struct edge {
    int u;          // G0 node
    int v;          // G0 node
    double key;        //
    double ckey;       // for soft heap
    int id;         // for soft heap delete
    bool removed;   // lazy removal from SoftHeap
    bool corrupted;

    friend bool operator<(edge e0, edge e1) {
        if (e0.key == e1.key) return e0.id < e1.id;
        return e0.key < e1.key;
    }

    std::string toString() const {
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
std::pair<int, std::vector<int>> tHierarchySize(int n, int m) {
    if (n < 3) {
        return {{3}, {10, 4, 1}};
    }
    // Arbitrary sequence for now that grows really fast -- no point in using Ackermann sequence.
    int desiredNumberOfLeaves[7] = {1, 1 << 2, 1 << 6, 1 << 16, 1 << 30, 1 << 30, 1 << 30};
    int desiredNumberOfChildren[7];
    desiredNumberOfChildren[0] = 1;
    for (int i = 1; i < 7; ++i) {
        desiredNumberOfChildren[i] = desiredNumberOfLeaves[i] / desiredNumberOfLeaves[i - 1];
    }
    std::vector<int> desiredSize;
    int d = 0;
    for (; desiredNumberOfLeaves[d] < n; ++d) {}
    for (int i = d; i >= 0; --i) desiredSize.push_back(desiredNumberOfChildren[i]);

    return {d, desiredSize};
}

NetworKit::Graph Chazelle2000MinimumSpanningTree::msf(NetworKit::Graph G, int t) {
    auto tree = NetworKit::GraphTools::copyNodes(G);
    auto connected_components = NetworKit::ConnectedComponents(G);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (auto i = 0; i < connected_components.numberOfComponents(); i++) {
        auto subgraph = NetworKit::GraphTools::subgraphFromNodes(G, [](const auto& components) {
            std::unordered_set<NetworKit::node> s;
            for (auto n : components) {
                s.insert(n);
            }
            return s;
        }(components[i]));

        using NetworKit::node;
        using NetworKit::edgeweight;
        std::unordered_map<node, node> mappingTo;
        std::unordered_map<node, node> mappingFrom;
        subgraph.forNodes([&](node u){
            int n = mappingTo.size();
            mappingTo[u] = n;
            mappingFrom[n] = u;
        });
        subgraph = NetworKit::GraphTools::getCompactedGraph(subgraph, mappingTo);
        auto tree1 = mst(subgraph, t);
        tree1.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
            tree.addEdge(mappingFrom[u], mappingFrom[v], ew);
        });
    }
    return tree;
}

namespace {
bool isConnected(NetworKit::Graph G) {
    auto connected_components = NetworKit::ConnectedComponents(G);
    connected_components.run();
    auto components = connected_components.getComponents();
    return components.size() == 1;
}
}

bool noMultiEdges(NetworKit::Graph G) {
    using NetworKit::node;
    std::set<std::pair<node, node>> s;
    bool ret = true;
    G.forEdges([&](node u, node v){
        if (s.contains(std::minmax(u, v))) ret = false;
        s.insert(std::minmax(u, v));
    });
    return ret;
}

const double LARGE_VALUE = __DBL_MAX__;

/**
 * Expects that:
 * G is connected
 * G is undirected
 */
NetworKit::Graph Chazelle2000MinimumSpanningTree::mst(NetworKit::Graph G, int t) {
    using std::minmax;
    using std::swap;
    using std::set;
    using std::pair;
    using std::map;
    using std::tuple;
    using std::vector;
    using NetworKit::Graph;
    using NetworKit::node;
    using NetworKit::edgeweight;

    assert(isConnected(G));
    assert(!G.isDirected());
    assert(noMultiEdges(G));

    // [STEP 1]
    if (t <= 1 || G.numberOfNodes() < MIN_NUMBER_NODES) {
        BoruvkaMinimumSpanningTree boruvka(G);
        boruvka.run();
        return boruvka.getForest();
    }

    // [STEP 2]
    auto [G0, edge_G0_to_G, forest] = boruvkaSteps(G, C);
    if (G0.numberOfEdges() == 0) {
        return forest;
    }

    edge dummy_edge{0, 0, __DBL_MAX__, __DBL_MAX__, 0, true, true};
    vector<vector<SoftHeap<edge*>>> heaps;
    set<NodePair> contractedEdges;
    set<NodePair> badEdges;
    vector<edge> edges(G0.numberOfEdges());
    vector<set<int>> Cz;
    vector<int> parent(G0.numberOfNodes());
    vector<bool> fusionNode(G0.numberOfNodes(), false);
    vector<vector<pair<NodePair, double>>> minLink;
    map<NodePair, int> G0_edge_id;
    vector<bool> visited(G0.numberOfNodes(), false);
    vector<std::stack<int>> border(G0.numberOfNodes(), std::stack<int>());
    assert(!G.isDirected());
    assert(!G0.isDirected());
    auto [d, desiredSize] = tHierarchySize(G0.numberOfNodes(), G0.numberOfEdges());


    // we have Cz0, Cz1, .... Czk
    auto k = [&](){return Cz.size() - 1;};

    auto makeEdge = [&](node u, node v, edgeweight ew) {
        return edge{u, v, ew, ew, G0_edge_id[minmax(u, v)], false, false};
    };

    auto minBorderEdge = [&](){
        double minKey = __DBL_MAX__;
        int mini = -1, minj = -1;

        for (int i = 0; i < heaps.size(); ++i) {
            for (int j = 0; j < heaps[i].size(); ++j) {
                auto ret = heaps[i][j].lookupMin();
                if (std::holds_alternative<bool>(ret)) continue;
                auto val = std::get<edge*>(ret);
                if (val->ckey < minKey) {
                    mini = i;
                    minj = j;
                    minKey = val->ckey;
                }
            }
        }
        edge* minEdge = heaps[mini][minj].extractMin();
        return std::tuple{*minEdge, 0, 0};
    };

    auto findAndDeleteEdgeFromHeaps = [&](node u, node v) {
        auto eid = G0_edge_id[minmax(u, v)];
        edges[eid].removed = true;
        if (edges[eid].corrupted) {
            badEdges.insert(minmax(u, v));
        }
        return edges[eid];
    };

    auto nodeTopParent = [&](node u) {
        int p = static_cast<int>(u);
        while (parent[p] != p) p = parent[p];
        return p;
    };

    auto CziOfNode = [&](node u) {
        int p = nodeTopParent(u);
        for (int i = 0; i < Cz.size(); ++i) {
            if (Cz[i].contains(p)) {
                return i;
            }
        }
        throw std::runtime_error("u node should be found");
    };

    auto insertNewBorderEdge = [&](edge* e, node u, node v) {
        int j = Cz.size();
        if (border[v].empty()) {
            heaps[0][j] = std::move(insert(std::move(heaps[0][j]), e));
        } else {
            int i = 1 + CziOfNode(border[v].top());
            heaps[i][j] = std::move(insert(std::move(heaps[i][j]), e));
        }
        border[v].push(u);
    };

    auto leftmostSmallerMinLink = [&](int ckey){
        NodePair ab;

        for (int i = 0; i < k(); ++i) {
            for (int j = i + 1; j <= k(); ++j) {
                if (minLink[i][j].second < 0) continue;
                if (minLink[i][j].second <= ckey) {
                    return tuple<int, int, NodePair>{i, j, minLink[i][j].first};
                }
            }
        }
        return tuple<int, int, NodePair>{-1, -1, ab};
    };

    set<edge> insideEdges;

    auto updateMinLinks = [&]() {
        int last_i = k();
        for (int i = 0; i < last_i; ++i) {
            int best_j = last_i;
            double best_val = minLink[i][last_i].second;
            for (int j = last_i + 1; j < minLink[i].size(); ++j) {
                if (minLink[i][j].second < 0) continue;
                if (minLink[i][j].second >= best_val && best_val >= 0) continue;
                best_j = j;
                best_val = minLink[i][j].second;
            }
            minLink[i][last_i] = minLink[i][best_j];
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                if (i >= k() || j > k()) {
                    minLink[i][j] = {{-1, -1}, __DBL_MAX__};
                }
            }
        }
    };

    auto popHeaps = [&]() {
        int k = Cz.size();
        heaps[k - 1][k - 1] = std::move(meld(std::move(heaps[k - 1][k - 1]),
            std::move(heaps[k - 1][k])));
        heaps[k - 1][k - 1] = std::move(meld(std::move(heaps[k - 1][k - 1]),
            std::move(heaps[k][k])));

        for (int i = 0; i < k - 1; ++i) {
            heaps[i][k - 1] = std::move(meld(std::move(heaps[i][k - 1]), std::move(heaps[i][k])));
        }
        heaps.pop_back();
        for (auto& h : heaps) {
            h.pop_back();
        }
    };

    auto retraction = [&]() {
        int newP = parent.size();
        parent.push_back(newP);
        fusionNode.push_back(false);
        auto vertices = Cz[Cz.size() - 1];
        popHeaps();
        Cz.pop_back();
        if (vertices.size() == 1) {
            parent.pop_back();
            fusionNode.pop_back();
            newP = *vertices.begin();
        } else {
            for (auto v : vertices) {
                parent[v] = newP;
            }
        }
        Cz[Cz.size() - 1].insert(newP);
        updateMinLinks();
    };

    auto fusion = [&]() -> edge {
        auto [uvEdge, mini, minj] = minBorderEdge();
        node uG0{uvEdge.u}, vG0{uvEdge.v};
        auto [link_i, link_j, ab] = leftmostSmallerMinLink(uvEdge.ckey);
        if (link_i == -1) return uvEdge;
        auto [a, b] = ab;
        auto ap = nodeTopParent(a);
        auto bp = nodeTopParent(b);
        if (!Cz[link_i].contains(ap)) {
            swap(a, b);
            swap(ap, bp);
        }

        int newP = parent.size();
        parent.push_back(newP);
        fusionNode.push_back(true);
        while (k() > link_i) {
            for (auto v : Cz[k()]) {
                parent[v] = newP;
            }
            popHeaps();
            Cz.pop_back();
        }
        Cz[k()].erase(ap);
        parent[ap] = newP;
        fusionNode[newP] = true;
        Cz[k()].insert(newP);

        updateMinLinks();
        return uvEdge;
    };

    auto updateMinLinksWithEdge = [&](edge e) {
        int i = CziOfNode(e.u), j = CziOfNode(e.v);
        if (i > j) swap(i, j);
        if (e.ckey < minLink[i][j].second) {
            minLink[i][j] = {minmax(e.u, e.v), e.ckey};
        }
    };

    int extensions = 0;
    auto extension = [&](edge e) {
        extensions += 1;
        node u = e.u, v = e.v;
        if (!visited[u]) swap(u, v);

        heaps.push_back(vector<SoftHeap<edge*>>());
        for (int i = 0; i <= Cz.size(); ++i) heaps[heaps.size() - 1].push_back(
            SoftHeap<edge*>(&dummy_edge, 0.1));
        Cz.push_back({v});
        for (int i = 0; i < heaps.size(); ++i) {
            heaps[i].push_back(SoftHeap<edge*>(&dummy_edge, 0.1));
        }
        visited[v] = true;
        std::vector<NodePair> newBorderEdges;
        G0.forEdgesOf(v, [&](node vv, node w, edgeweight ew) {
            if (visited[w]) {
                edge e = findAndDeleteEdgeFromHeaps(v, w);
                updateMinLinksWithEdge(e);
                insideEdges.insert(e);
            } else {
                int eid = G0_edge_id[minmax({v, w})];
                insertNewBorderEdge(&edges[eid], v, w);
            }
        });
        updateMinLinks();
    };

    auto shouldRetract = [&]() {
        return Cz.size() == d || (Cz[k()].size() >= desiredSize[k()] && k() > 0);
    };

    auto shouldFinish = [&]() {
        return extensions >= G0.numberOfNodes() - 1;
    };

    auto initialization = [&]() {
        int eidNow = 0;
        G0.forEdges([&](node u, node v, edgeweight ew) {
            edges[eidNow] = edge{u, v, ew, ew, eidNow, false, false};
            G0_edge_id[minmax(u, v)] = eidNow;
            eidNow += 1;
        });

        for (int i = 0; i < parent.size(); i++) parent[i] = i;
        Cz.push_back({0});
        heaps.push_back({SoftHeap<edge*>(&dummy_edge, 0.1), SoftHeap<edge*>(&dummy_edge, 0.1)});
        heaps.push_back({SoftHeap<edge*>(&dummy_edge, 0.1), SoftHeap<edge*>(&dummy_edge, 0.1)});
        G0.forEdgesOf(0, [&](node u, node v, edgeweight ew) {
            int eid = G0_edge_id[minmax(u, v)];
            insertNewBorderEdge(&edges[eid], u, v);
        });
        minLink = vector<vector<pair<NodePair, double>>>(d,
            vector<pair<NodePair, double>>(d, {{-1, -1}, __DBL_MAX__}));
        visited[0] = true;
    };

    initialization();
    while (!shouldFinish()) {
        if (shouldRetract()) {
            retraction();
        } else {
            edge e = fusion();
            extension(e);
        }
    }

    while (Cz.size() > 1) {
        int newP = parent.size();
        parent.push_back(newP);
        fusionNode.push_back(false);
        auto vertices = Cz[Cz.size() - 1];

        heaps[k() - 1][k() - 1] = std::move(
            meld(std::move(heaps[k()][k()]), std::move(heaps[k() - 1][k() - 1])));
        heaps.pop_back();
        Cz.pop_back();

        for (auto v : vertices) {
            parent[v] = newP;
        }
        Cz[Cz.size() - 1].insert(newP);
    }
    int root = parent.size();
    parent.push_back(root);
    fusionNode.push_back(false);
    for (int v : Cz[0]) {
        parent[v] = root;
    }

    for (int i = 0; i < edges.size(); ++i) {
        if (!edges[i].corrupted) continue;
        badEdges.insert(minmax(edges[i].u, edges[i].v));
    }

    // [STEP 4]
    vector<map<node, node>> fromG0Maps(parent.size());
    vector<map<NodePair, NodePair>> toG0edgeMaps(parent.size());
    vector<int> depth(parent.size(), 0);

    for (int i = 0; i < G0.numberOfNodes(); ++i) {
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
            int new_v = fromG0Maps[v].size();
            fromG0Maps[v].insert({prev_v, new_v});
        }
    }

    vector<Graph> CzGraphs(parent.size(), Graph(0, true));
    for (int i = G0.numberOfNodes(); i < parent.size(); ++i) {
        CzGraphs[i] = Graph(fromG0Maps[i].size(), true);
    }

    auto lcaLR = [&](int l, int r) {
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

        return tuple<int, int, int>{l, prev_l, prev_r};
    };

    vector<map<NodePair, double>> smallestEdge(CzGraphs.size());

    G0.forEdges([&](node u, node v, edgeweight ew) {
        if (badEdges.contains(minmax(u, v))) return;
        auto [lca, l, r] = lcaLR(u, v);
        auto minorL = fromG0Maps[lca][l];
        auto minorR = fromG0Maps[lca][r];
        auto np = minmax(minorL, minorR);
        if (!smallestEdge[lca].contains(np)) {
            smallestEdge[lca][np] = ew;
            toG0edgeMaps[lca][np] = minmax(u, v);
        }
        if (smallestEdge[lca][np] > ew) {
            smallestEdge[lca][np] = ew;
            toG0edgeMaps[lca][np] = minmax(u, v);
        }
    });
    for (int i = G0.numberOfNodes(); i < CzGraphs.size(); ++i) {
        for (auto [key, val] : smallestEdge[i]) {
            CzGraphs[i].addEdge(key.first, key.second, val);
        }
    }

    Graph G0copy = G0;
    BoruvkaMinimumSpanningTree bbb(G0copy);
    bbb.run();
    auto actualMst = bbb.getForest();
    set<NodePair> actualMstEdges;
    actualMst.forEdges([&](node u, node v, edgeweight ew) {
        actualMstEdges.insert(minmax(u, v));
    });

    for (int i = G0.numberOfNodes(); i < CzGraphs.size(); ++i) {
        if (fusionNode[i]) continue;
        CzGraphs[i] = msf(CzGraphs[i], t - 1);
    }

    auto F = NetworKit::GraphTools::copyNodes(G0);
    for (int i = G0.numberOfNodes(); i < CzGraphs.size(); ++i) {
        CzGraphs[i].forEdges([&](node u, node v){
            auto [x, y] = toG0edgeMaps[i][minmax(u, v)];
            F.addEdge(x, y, edges[G0_edge_id[minmax(x, y)]].key);
        });
    }

    for (auto [x, y] : badEdges) {
        F.addEdge(x, y, edges[G0_edge_id[minmax(x, y)]].key);
    }

    // [STEP 5]
    Graph res = mst(F, t);

    res.forEdges([&](node u, node v, edgeweight ew) {
        auto [uG, vG] = edge_G0_to_G[minmax(u, v)];
        forest.addEdge(uG, vG, ew);
    });

    for (int i = 0; i < parent.size(); ++i) {
        int v = i;
        while (v != parent[v]) {
            v = parent[v];
        }
    }

    return forest;
}

}  /* namespace Koala */
