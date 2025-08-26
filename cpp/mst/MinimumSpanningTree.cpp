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
        // TODO It would be nice to clear the E map of the edges that actually don't exist in here
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
    std::uniform_int_distribution<int> distrib(0, graph->numberOfNodes() - 1);

    for (int i = 0; i < C / eps; ++i) {
        int v = distrib(generator);
        d_est = std::max(d_est, static_cast<int>(graph->degreeOut(v)));
    }
    return d_est;
}

float ChazelleRubinfeldTrevisanMinimumSpanningTree::calculateApproximateCCsCount(float eps, int bfs_bound, unsigned int w, unsigned int w_bound) const {
    int n = graph->numberOfNodes();
    int r = 1/eps/eps + 1;
    float estimate = 0;

    float d_est = calculateApproximateDegree(eps);
    std::uniform_int_distribution<int> distrib(0, n-1);

    for (int i = 0; i < r; ++i) {
        int u = distrib(generator);
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

        graph->forEdgesOf(u, [this, w_bound, &d_actual](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid id) {
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
        graph->forEdgesOf(u, [this, w_bound, &visited_edges, &q](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid id){
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

                graph->forEdgesOf(v, [this, w_bound, &visited_edges, &q, &vis](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid id) {
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
                    b = powl(2, coin_flips) * static_cast<float>(d_actual) / static_cast<float>(visited_edges);
                }
                break;
            }
        }

        estimate += b;
    }

    return estimate * static_cast<float>(n) / 2.0f / r;
}

float ChazelleRubinfeldTrevisanMinimumSpanningTree::calculateApproximateTreeWeight(float eps, unsigned int w) const {
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
    throw std::runtime_error("getForest method not supported for approximate MST algorithm. Try getTreeWeight");
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
}

void Chazelle2000MinimumSpanningTree::run() {
    // TODO Find proper starting t value
    // TODO Connected components in here
    auto G = graph.value();
    std::cout << "RUNNNING MST CHAZELLE ALGORITHM" << std::endl;
    // auto connected_components = NetworKit::ConnectedComponents(G);
    // connected_components.run();
    // auto components = connected_components.getComponents();
    // for (auto i = 0; i < connected_components.numberOfComponents(); i++) {
    //     auto subgraph = NetworKit::GraphTools::subgraphFromNodes(G, [](const auto& components){
    //         std::unordered_set<NetworKit::node> s;
    //         for (auto n: components){
    //             s.insert(n);
    //         }
    //         return s;
    //     }(components[i]));
    //     auto tree1 = mst(subgraph, 10);
    //     tree1.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
    //         tree->addEdge(u, v, ew);
    //     });
    // }
    
    tree = mst(G, 10);
    std::cout << "RUNNNING MST CHAZELLE ALGORITHM" << std::endl;
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
std::tuple<NetworKit::Graph, std::map<Koala::MinimumSpanningTree::NodePair, Koala::MinimumSpanningTree::NodePair>, NetworKit::Graph>
    Chazelle2000MinimumSpanningTree::boruvkaSteps(NetworKit::Graph G, int c) {
    assert(c >= 1);

    std::cout << "INSIDE BORVUKA STEPS" << std::endl;
    auto F = NetworKit::GraphTools::copyNodes(G);
    auto uf = NetworKit::UnionFind(graph->upperNodeIdBound());
    std::map<NodePair, NodePair> gMinor_to_G;
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        gMinor_to_G.insert({std::minmax(u, v), {u, v}});
    });
    auto gMinor = G;
    BoruvkaMinimumSpanningTree::iterate(gMinor, F, uf, gMinor_to_G, c, false);
    
    using NetworKit::node; using NetworKit::edgeweight;
    // gMinor.forEdges([&](node u, node v, edgeweight ew){
    //     auto [uu, vv] = gMinor_to_G[std::minmax(u,v)];
    //     std :: cout << "gMinor before: " << u << " " << v << " " << ew <<  ", " << uu << " " << vv << std::endl;
    // });
    auto map1 = NetworKit::GraphTools::getContinuousNodeIds(gMinor);
    auto gMinorCompact = NetworKit::GraphTools::getCompactedGraph(gMinor, map1);
    std::map<NodePair, NodePair> gMinorRemapped_to_G;
    // for (auto [key, val]: map1) {
    //     std :: cout << key << " -> " << val << std::endl; 
    // }
    // for (auto [key, val]: gMinor_to_G) {
    //     gMinorRemapped_to_G[std::minmax(map1[key.first], map1[key.second])] = val;
    // }
    for (auto const& [key, val] : gMinor_to_G) {
        auto oldU = key.first;
        auto oldV = key.second;
        auto newU_it = map1.find(oldU);
        auto newV_it = map1.find(oldV);

        if (newU_it == map1.end() || newV_it == map1.end()) {
            // std::cout << "[WARNING] Node not found in map1: "
            //         << oldU << " or " << oldV << std::endl;
            continue;
        }

        auto newU = newU_it->second;
        auto newV = newV_it->second;

        gMinorRemapped_to_G[std::minmax(newU, newV)] = val;

        // std::cout << "Mapping gMinor_to_G edge: (" 
        //         << oldU << ", " << oldV << ")  -> compacted ("
        //         << newU << ", " << newV << ")"
        //         << "   original G edge: (" << val.first << ", " << val.second << ")"
        //         << std::endl;
    }
    // gMinorCompact.forEdges([&](node u, node v, edgeweight ew){
    //     auto [uu, vv] = gMinorRemapped_to_G[std::minmax(u,v)];
    //     std :: cout << "gMinorCompact: " << u << " " << v << " " << ew <<  ", " << uu << " " << vv << std::endl;
    // });
    return {gMinorCompact, gMinorRemapped_to_G, F};
    
    // std::map<NetworKit::node, NetworKit::node> map1;
    // int nodeId = 0;
    // gMinor.forNodes([&](NetworKit::node u){
    //     map1[u] = nodeId;
    //     nodeId += 1;
    // });

    // auto gMinorRemapped = NetworKit::GraphTools::getRemappedGraph(gMinor, nodeId, [&](NetworKit::node u){return map1[u];});
    // for (auto [key, val]: gMinor_to_G) {
    //     gMinorRemapped_to_G[std::minmax(map1[key.first], map1[key.second])] = val;
    // }

    // return {gMinorRemapped, gMinorRemapped_to_G, F};
}

namespace{
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
}

// NetworKit::Graph Chazelle2000MinimumSpanningTree::msf2(NetworKit::Graph G, std::map<NodePair, int>& edgeId, int t) {
//     using std::minmax;
    
//     std::cout << "HUB" << std::endl;
//     // auto ID = G.getEdgeIntAttribute("ID");
//     // std::cout << "HUB 2" << std::endl;
//     // std::cout << ID[0] << std::endl;

//     // [STEP 1]
//     if( t <= 1 || G.numberOfNodes() < MIN_NUMBER_NODES) {
//         BoruvkaMinimumSpanningTree boruvka(G);
//         boruvka.run();
//         return boruvka.getForest();
//     }

    
//     // [STEP 2]
//     auto [G0, edgeId0, forest] = boruvkaSteps(G, edgeId, C);
//     if (G0.numberOfEdges() == 0) {
//         return forest;
//     }
    
//     // [STEP 3]
//     // [TODO]

//     /**
//      * I actually don't use the minors as I would also need to store
//      * a map for each minor that maps the edges, will be easier
//      * to calculate MSF recursively on the fly and add edges to F
//      * Already commited Cz[i] graph representing a minor of a subgraph of G0
//      */
//     // std::vector<std::pair<NetworKit::Graph, std::map<NodePair, int>>> minors;

//     /**
//      * Using binomial heaps for now as they have the `update` method
//      *  I'll also need a meld method but will implement it naively for now
//      * 
//      * [TODO] Change it to SoftHeap
//      */
//     std::vector<std::vector<Koala::BinomialHeap<edge>>> heaps;
    
//     /**
//      * Cz[i] holds the vertices of GT that are in Cz[i] minor subgraph
//      * maybe some more information will be needed, let's keep it simple for now
//      */
//     std::vector<std::set<NetworKit::node>> Cz;
//     // std::vector<std::pair<NetworKit::Graph, std::map<NodePair, int>>> Cz;
    
//     // [TODO] figure out where will minLink and chainLink be used...
//     /**
//      * min link should be a GT edge
//      */
//     std::vector<std::vector<edge>> minLink;
//     // std::vector< edgeId > chain_link;
//     // std::vector<std::vector<std::vector<edge>>> CzEdges;

//     /**
//      * In F we store the edges contracted during the [STEP 3] 
//      * They may not end up in the actual MSF returned by this function,
//      * as the [STEP 5] calculates MSF based on graphs F and B
//      */
//     NetworKit::Graph F = NetworKit::GraphTools::copyNodes(G0);
//     /**
//      * Graph of BAD edges created during [STEP 3]
//      * Need to be reprocessed in [STEP 5] to get a final MSF
//      */
//     NetworKit::Graph B = NetworKit::GraphTools::copyNodes(G0);
//     /**
//      * GT represents a graph of the current active path in T hierarchy
//      * i.e. each Cz has its contracted vertices which are connected
//      * by some original edges from the G0 graph
//      */

//     NetworKit::Graph GT;
//     std::map<NodePair, int> GT_eid;
    
//     /**
//      * maps unique edge id to (u, v) edge in G0 graph
//      */
//     std::map<int, NodePair> eid_2_G0;
//     G0.forEdges([&](auto u, auto v){
//         eid_2_G0[edgeId0[minmax(u, v)]] = minmax(u, v);
//     });

//     std::map<NetworKit::node, int> G0_GT;
//     for (auto v : G0.nodeRange()) {
//         G0_GT[v] = -1;
//     }

//     /**
//      * Performs a retraction operation defined by Chazelle
//      * force flag indicates that we want to retract the only left Cz
//      * i.e. when the Cz.size() == 1. Used to assert proper algo behaviour
//      */
//     auto retraction = [&](bool force = false){
//         // We need to explicitly call retraction with force 
//         // if we want to retract the last Cz component
//         assert(!force || Cz.size() > 1);
        
//         const auto& V =  Cz[Cz.size() - 1];
//         std::vector<std::pair<NodePair, NetworKit::edgeweight>> E;
//         NetworKit::Graph GMinor;
//         std::map<NetworKit::node, NetworKit::node> V_2_GMinorV;
//         std::map<NodePair, int> GMinorEdgeId;
//         std::map<NetworKit::node, std::pair<NodePair, NetworKit::edgeweight>> smallestEdge;

//         for (NetworKit::index v_id : V) {
//             V_2_GMinorV[v_id] = GMinor.addNode();
//             GT.forEdgesOf(v_id, [&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
//                 E.push_back({{u, v}, ew});
//             });
//         }

//         for (auto[e, ew] : E) {
//             auto[u, v] = e;

//             // uv edge is in the GMinor 
//             if (V.contains(u) && V.contains(v)) {
//                 auto uGM = V_2_GMinorV[u];
//                 auto vGM = V_2_GMinorV[v];
//                 // [THINK TODO] is ew correct here?
//                 GMinor.addEdge(uGM, vGM, ew);
//                 // [TODO] Check that edgeId0 is the map that maps u, v from GT (same as G0) to proper edgeId in G
//                 GMinorEdgeId[std::minmax(uGM, vGM)] = edgeId0[std::minmax(u, v)];
//             } else {
//                 if (!V.contains(v)) {
//                     std::swap(v, u);
//                 }
//                 // v in V
//                 // u in other Czs
//                 if (!smallestEdge.contains(u)) {
//                     smallestEdge[u] = {{u, v}, ew};
//                 }
//                 if (ew < smallestEdge[u].second) {
//                     smallestEdge[u] = {{u, v}, ew};
//                 }
//             }
//         }

//         // All neccessary information was extracted from the GT
//         // Now need to transform the GT itself
//         // [THINK TODO] What about the border structure ???

//         Cz.pop_back();
//         if (Cz.size() == 0) {
//             // [THINK TODO] That's only happening when contracting the last minor of the whole graph
//             return;
//         }
//         for (auto v : V) {
//             GT.removeNode(v);
//         }

//         NetworKit::node newNode = GT.addNode();
//         Cz[Cz.size() - 1].insert(newNode);
//         for (auto [u, e] : smallestEdge) {
//             GT.addEdge(newNode, u, e.second);
//         }
        
//         // GT and Cz are now in a valid state

//         // [STEP 4]
//         // Calculate the msf of GMinor now, why not
//         // [TODO] Need to remove the bad edges in here or during the construction of GMinor
//         // [TODO] Check whether the ew is correct or incorrect...
//         auto minorForest = msf(GMinor, GMinorEdgeId, t - 1);
//         minorForest.forEdges([&](auto u, auto v, NetworKit::edgeweight ew){
//             int eid = GMinorEdgeId[std::minmax(u, v)];
//             // auto [uF, vF] = eid_2_G0[eid];
//             // F.addEdge(uF, vF, ew);
//         });

//         /**
//          * [THINK] After contraction of the last Cz to a vertex there will be multiple
//          * edges in the GT graph, figure out what to do with that, it may not be obvious
//          * to leave just the smallest one.
//          * 
//          * Leaving the smallest one for now
//          */

//          /**
//           * [THINK TODO] BORDER STRUCTURE will be a mess
//           * Maybe a map from G0 vertices to GT vertices will be needed ???
//           */
//     };

//     auto minBorderEdge = [&]() {
//         int minVal = INT32_MAX;
//         int mini = -1, minj = -1;
//         for (int i = 0; i < heaps.size(); ++i) {
//             for (int j = 0; j < heaps[i].size(); ++j) {
//                 if ( heaps[i][j].top().val < minVal ) {
//                     mini = i;
//                     minj = j;
//                     minVal = heaps[i][j].top().val;
//                 }
//             }
//         }

//         edge minEdge = heaps[mini][minj].top();
//         return std::tuple{minEdge, mini, minj};
//     };

//     auto leftmostSmallerMinLink = [&](int val) {
//         for (int i = 0; i < minLink.size(); ++i) {
//             for (int j = 0; j < minLink[i].size(); ++j) {
//                 if( minLink[i][j].val <= val ) {
//                     return std::tuple{true, i, j};
//                 }
//             }
//         }
//         return std::tuple{false, -1, -1};
//     };

//     auto fusion = [&](int linki, int linkj) {
//         auto e = minLink[linki][linkj];
//         NetworKit::node a = e.u, b = e.v;
//         if (!Cz[linki].contains(a)) {
//             a = e.v;
//             b = e.u;
//         }
//         assert(Cz[linki].contains(a) && Cz[linkj].contains(b));

//         std::set<NetworKit::node> contractedV;
//         std::vector<std::pair<NodePair, NetworKit::edgeweight>> E;
//         std::map<NetworKit::node, std::pair<NodePair, NetworKit::edgeweight>> smallestEdge;

//         contractedV.insert(a);
//         for (int i = linki + 1; i < Cz.size(); ++i) {
//             for (auto v : Cz[i]) {
//                 contractedV.insert(v);
//                 GT.forEdgesOf(v, [&](auto u, auto v, NetworKit::edgeweight ew){
//                     E.push_back({minmax(u, v), ew});
//                 });
//             }
//         }

//         GT.forEdgesOf(a, [&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
//             if (contractedV.contains(u) && contractedV.contains(v)) {
//                 return;
//             }
//             if(u != a) {
//                 std::swap(u, v);
//             }
//             smallestEdge[v] = {minmax(u, v), ew};
//         });

//         for (auto e : E) {
//             auto [uv, ew] = e;
//             auto [u, v] = uv;

//             // Need to contract that edge
//             if (contractedV.contains(u) && contractedV.contains(v)) {
//                 auto[uF, vF] = eid_2_G0[GT_eid[minmax(u, v)]];
//                 // [TODO] what ew is good in here?
//                 F.addEdge(uF, vF, ew);
//             } else {
//                 if (!contractedV.contains(u)) {
//                     std::swap(u, v);
//                 }
                
//                 if (ew < smallestEdge[v].second) {
//                     smallestEdge[v] = {minmax(u, v), ew};
//                 }
//             }
//         }

//         for (auto v : contractedV) {
//             if (v == a) {
//                 continue;
//             }

//             GT.removeNode(v);
//         }

//         // Remove old edges from a
//         GT.removeAdjacentEdges(a, [](auto x){return true;}, false);
//         GT.removeAdjacentEdges(a, [](auto x){return true;}, true);
//         // Add new edges to a
//         // [TODO] need to update edge maps in here
//         for (auto [v, e] : smallestEdge) {
//             GT.addEdge(a, v, e.second);
//         }

//         // [TODO] deal with the heaps...
//     };

//     auto extension = [&]() {
//         {
//             auto [minEdge, i, j] = minBorderEdge();
//             auto [needFusion, linki, linkj] = leftmostSmallerMinLink(minEdge.val);
    
//             if ( needFusion ) {
//                 fusion(linki, linkj);
//             }
//         }
//         auto [minEdge, i, j] = minBorderEdge();
//         heaps[i][j].pop();

//         NetworKit::node bGT = GT.addNode();
//         Cz.push_back({bGT});
        
//         // auto bG0 = G0_GT.conatins()
//         // G0.forEdgesOf(bG0, [&](auto u, auto v, NetworKit::edgeweight ew){
//         //     if (u != b) std::swap(u, v);
//         //     if(G0_GT.contains(v)) {
                
//         //     } 
//         // });
//         // [TODO] add proper edges to GT
//         // [TODO] add minLinks
//         // [TODO] chain link == minEdge ???
//         // [TODO] add new edges to heaps
//     };

//     // [TODO]
//     auto heapsEmpty = [&]() {
//         return false;
//     };

//     // [TODO]
//     auto shouldRetract = [&]() {
//         return false;
//     };

//     while ( !heapsEmpty() ) {
//           if ( shouldRetract() ) {
//               retraction();
//           } else {
//               extension();
//          }
//      }
     

//     // [STEP 4]
    
//     /**
//      * Call recursivly msf on each minor from minors
//      * let F be the sum of all of resulting forests
//      * As in [STEP 5] we will call msf on F and on Bad edges,
//      * its important not to put these edges to the actual msf. 
//      * 
//      * It's actually done during the retraction...
//      */

//     // NetworKit::Graph F = NetworKit::GraphTools::copyNodes(G0);
//     // for (auto minor : minors) {
//     //     /**
//     //      * [TODO]
//     //      * F0 = msf(minor \ B, t - 1);
//     //      * [TODO]
//     //      * F += F0
//     //      */
//     // }
    

//     // [STEP 5]
//     // [TODO]
//     // return msf( F + B, t);

// }
namespace{
    int pow2Clamped(int b) {
        if (b > 30) {
            return 1 << 30;
        } else {
            return 1 << b;
        }
    }
}

int Chazelle2000MinimumSpanningTree::verticesOnLevel(int t) {
    // it should be some kind of ackermann function
    // but for now let's make it 2**(2**n)
    // it's easier to compute and actually has some small values for the first couple of values.
    return pow2Clamped(pow2Clamped(t));
}

// /**
//  * Expects that:
//  * G is connected
//  * G is undirected
//  */
// NetworKit::Graph Chazelle2000MinimumSpanningTree::mst(NetworKit::Graph G, int t) {
//     using std::minmax;
//     using std::swap;
//     using std::set;
//     using std::pair;
//     using std::map;
//     using std::tuple;
//     using std::vector;
//     using NetworKit::Graph;
//     using NetworKit::node;
//     using NetworKit::edgeweight;
    
//     std::cout << "INSIDE CHAZELLE" << std::endl;
//     // [STEP 1]
//     if (t <= 1 || G.numberOfNodes() < MIN_NUMBER_NODES) {
//         BoruvkaMinimumSpanningTree boruvka(G);
//         boruvka.run();
//         return boruvka.getForest();
//     }
    
//     // [STEP 2]
//     std::cout << "BEFORE BORUVKA STEPS" << std::endl;
//     auto [G0, edge_G0_to_G, forest] = boruvkaSteps(G, C);
//     if (G0.numberOfEdges() == 0) {
//         std::cout << "BORUVKA STEPS FINISHED" << std::endl;
//         // assert(false);
//         return forest;
//     }

//     // [STEP 3]
//     /* G0 nodes */
//     // vector<vector<Koala::BinomialHeap<edge*>>> heaps;
//     // for now let's make it work with only 1 heap
//     edge dummy_edge;
//     SoftHeap<edge*> heaps(&dummy_edge, 0.1);
//     heaps.extractMin();
//     Graph F = NetworKit::GraphTools::copyNodes(G0);
//     Graph B = NetworKit::GraphTools::copyNodes(G0);
//     vector<edge> edges(G0.numberOfEdges());
//     /* GT nodes */
//     Graph GT;
//     vector<set<node>> Cz;
//     vector<vector<pair<NodePair, double>>> minLink;
//     // vector<vector<double>> minLink;

//     vector<NodePair> chainLink;
//     /* maps */
//     map<NodePair, NodePair> edge_GT_to_G0;
//     // map<int, NodePair> eid_2_G0;
//     map<node, node> node_G0_to_GT;
//     map<NodePair, int> G0_edge_id;

//     assert(!G.isDirected());
//     assert(!G0.isDirected());
//     assert(!GT.isDirected());

//     int eidNow = 0;
//     G0.forEdges([&](node u, node v, edgeweight ew){
//         edges[eidNow] = edge{u, v, ew, ew, eidNow, false};
//         G0_edge_id[minmax(u, v)] = eidNow;
//         eidNow += 1;
//     });

//     auto removeGTNode = [&](node u) {
//         GT.forEdgesOf(u, [&](node uGT, node vGT){
//             edge_GT_to_G0.erase(minmax(uGT, vGT));
//         });
//         GT.removeNode(u);
//     };

//     // Returns a tuple {minEdge, min_i, min_j} so that heaps[min_i][min_j].top() == minEdge;
//     auto minBorderEdge = [&](){
//         int minKey = INT32_MAX;
//         int mini = -1, minj = -1;
//         // [TODO heaps]
//         // for (int i = 0; i < heaps.size(); ++i) {
//         //     for (int j = 0; j < heaps[i].size(); ++j) {
//         //         // lazy delete, will work well with soft heaps
//         //         while (!heaps[i][j].empty() && heaps[i][j].top()->removed) {
//         //             heaps[i][j].pop();
//         //         }
//         //         if (!heaps[i][j].empty() && heaps[i][j].top()->key < minKey ) {
//         //             mini = i;
//         //             minj = j;
//         //             minKey = heaps[i][j].top()->key;
//         //         }
//         //     }
//         // }
//         // edge minEdge = *heaps[mini][minj].top();
//         // return std::tuple{minEdge, mini, minj};
//         edge* minEdge = heaps.extractMin();
//         while(minEdge->removed) {
//             minEdge = heaps.extractMin();
//         }
//         return std::tuple{minEdge, 0, 0};
//     };

//     auto findAndDeleteEdgeFromHeaps = [&](NodePair edgeG0) {
//         // [TODO] Implement
//         auto eid = G0_edge_id[minmax(edgeG0)];
//         edges[eid].removed = true;
//         return edges[eid];
//     };

//     auto insertNewBorderEdges = [&](const vector<edge>& newBorderEdges) {
//         // [TODO] Implement
//         for (auto e: newBorderEdges) {
//             NodePair np = minmax{e.u, e.v};
//             auto eid = G0_edge_id[np];
//             edges[eid] = e;
//             heaps = std::move(insert(std::move(heaps), &edges[eid]));
//         }
//         return;
//     };

//     // [TODO] Check all edge{} constructions and change to this...
//     auto makeEdge = [&](node u, node v, edgeweight ew) {
//         return edge{u, v, ew, ew, G0_edge_id[minmax(u, v)], false};
//     };

//     auto leftmostSmallerMinLink = [&](int ckey){
//         int i = 0;
//         int j = 0;
//         NodePair abGT;

//         for (; i < minLink.size(); ++i) {
//             for (; j < minLink[i].size(); ++j) {
//                 // doesnt exist
//                 if (minLink[i][j].second < 0) continue;
//                 if (minLink[i][j].second <= ckey) {
//                     return tuple<int, int, NodePair>{i, j, minLink[i][j].first};
//                 }
//             }
//         }

//         return tuple<int, int, NodePair>{-1, -1, abGT};
//     };

//     auto contractGTedge = [&](node u, node v){
//         auto [uG0, vG0] = edge_GT_to_G0[minmax(u, v)];
//         int key = edges[G0_edge_id[minmax(uG0, vG0)]].key;
//         F.addEdge(uG0, vG0, key);
//     };

//     auto updateMinLinks = [&](int last_i) {
//         for (int i = 0; i < last_i; ++i) {
//             int best_j = last_i;
//             double best_val = minLink[i][last_i].second;
//             for (int j = last_i + 1; j < minLink[i].size(); ++j) {
//                 if (minLink[i][j].second < 0) continue;
//                 if (minLink[i][j].second >= best_val) continue;
//                 best_j = j;
//                 best_val = minLink[i][j].second;
//             }
//             minLink[i][last_i] = minLink[i][best_j];
//         }
//         minLink.resize(last_i + 1);
//     };

//     auto updateHeaps = [&](int last_i) {
//         // [TODO] when heaps
//     };

//     auto fusion = [&](){
//         auto [uvEdge, mini, minj] = minBorderEdge();
//         node uG0{uvEdge->u}, vG0{uvEdge->v};
//         auto [link_i, link_j, abGT] = leftmostSmallerMinLink(uvEdge->ckey); // [TODO CHECK] ckey ???
//         if (link_i == -1) return uvEdge;
        
        
//         auto [aGT, bGT] = abGT;
//         // [TODO] check wheter we need to swap a and b, this may apply to every edge retrieval...

//         set<node> nodesToFuse;
//         nodesToFuse.insert(aGT);
//         // Better insert it twice than never
//         nodesToFuse.insert(bGT);
//         vector<pair<NodePair, edgeweight>> edgesToProcess;
//         map<node, edge> smallestEdge;
//         for (int i = link_i + 1; i < Cz.size(); ++i) {
//             for (node w : Cz[i]) {
//                 nodesToFuse.insert(w);
//                 GT.forEdgesOf(w, [&](node w, node x, edgeweight ew){
//                     edgesToProcess.push_back({{w, x}, ew});
//                 });
//             }
//         }
//         GT.forEdgesOf(aGT, [&](node aGT2, node xGT, edgeweight ew) {
//             assert(aGT == aGT2);
//             if (nodesToFuse.contains(xGT)) return;
//             smallestEdge[xGT] = edge {
//                 aGT, xGT, ew, ew
//             };
//         });

//         for (auto [wxGT, ew] : edgesToProcess) {
//             auto [wGT, xGT] = wxGT;
//             if (nodesToFuse.contains(xGT)) swap(wGT, xGT);
//             if (nodesToFuse.contains(xGT)) {
//                 // contract edge
//                 contractGTedge(wGT, xGT);
//             } else {
//                 if (!smallestEdge.contains(xGT)) {
//                     smallestEdge[xGT] = edge {
//                         wGT, xGT, ew, ew
//                     };
//                 }

//                 if (ew < smallestEdge[xGT].key) {
//                     smallestEdge[xGT] = edge {
//                         wGT, xGT, ew, ew
//                     };
//                 }
//             }
//         }

//         GT.removeAdjacentEdges(aGT, [](node _){return true;}, false);
//         GT.removeAdjacentEdges(aGT, [](node _){return true;}, true);
//         for (auto [xGT, axGT]: smallestEdge) {
//             GT.addEdge(aGT, xGT, axGT.key);
//             // [TODO] check if we need to add xa edge also...
//         }

//         updateMinLinks(link_i);

//         // [TODO heaps] 
//         // updateHeaps(link_i);

//         // [TODO CHECK] of by one error ??? And what about the complexity ???
//         chainLink.resize(link_i); 
//         for (node uGT : nodesToFuse) {
//             if (uGT == aGT) {
//                 continue;
//             }
//             removeGTNode(uGT);
//             // [TODO] check wheter removal works correctly
//             // GT.removeNode(uGT);
//         }
//         // [TODO CHECK] Same as for chain links
//         Cz.resize(link_i + 1);
//         return uvEdge;
//     };

//     auto extension = [&](){
//         // [TODO] assert(uv.u inside GT && uv.v outside GT)
//         auto [uvEdge, mini, minj] = minBorderEdge();
//         node uG0{uvEdge->u}, vG0{uvEdge->v};
//         node uGT{node_G0_to_GT[uG0]};

//         node vGT = GT.addNode();
//         Cz.push_back(set{vGT});
//         minLink.push_back(vector<pair<NodePair, double>>());
//         for (auto& v: minLink) {
//             v.push_back({{-1, -1}, -1});
//         }
//         // [TODO] update minLinks - resize
//         // [TODO] check whether that's correct for chainlink
//         // chainLink.push_back(makeEdge(uG0, vG0, uvEdge->ckey));
//         node_G0_to_GT[vG0] = vGT;

//         map<node, edge> smallestIncidentEdge;
//         std::vector<edge> newBorderEdges;
//         G0.forEdgesOf(vG0, [&](node vG0, node wG0, edgeweight ew){
//             if (node_G0_to_GT.contains(wG0)) {
//                 node wGT = node_G0_to_GT[wG0];
//                 // border edge, need to delete it from the heap it's in
//                 // edges[G0_edge_id[minmax(v, w)]].removed = true;
//                 edge e = findAndDeleteEdgeFromHeaps(minmax(vG0, wG0));
//                 if (!smallestIncidentEdge.contains(wGT)) {
//                     smallestIncidentEdge[wGT] = e;
//                 } else if (e.ckey < smallestIncidentEdge[wGT].ckey) {
//                     // [TODO] compare ckey or key?
//                     smallestIncidentEdge[wGT] = e;
//                 } else {
//                     // [TODO] can we just ignore that edge ???
//                 }
//             } else {
//                 // w lays outside GT => need to add it to a heap
//                 newBorderEdges.push_back(makeEdge(vG0, wG0, ew));
//             }
//         });
        
//         insertNewBorderEdges(newBorderEdges);
//         for (auto [wGT, e] : smallestIncidentEdge) {
//             // [TODO CHECK] In GT we can store edge with OG key
//             // as it will only be used in the recursion
//             GT.addEdge(vGT, wGT, e.key);
            
//             // [TODO] update minLinks - e is a min link
//             int wCzi = -1;
//             for (int i = 0; i < Cz.size(); ++i) {
//                 if (Cz[i].contains(wCzi)) {
//                     wCzi = i;
//                     break;
//                 }
//             }
//             assert(wCzi != -1);
//             if (minLink[wCzi][Cz.size() - 1] > e.ckey) {
//                 minLink[wCzi][Cz.size() - 1] = e.ckey;
//             }
//             // [TODO] maybe we need to also add wv edge?
//         }
//     };

//     auto minorMst = [&](const Graph& gMinor, const map<NodePair, NodePair>& edge_GMinor_to_G0){
//         Graph minorTree = mst(gMinor, t - 1);
//         minorTree.forEdges([&](node uMinor, node vMinor, edgeweight ew) {
//             auto [uG0, vG0] = edge_GMinor_to_G0[minmax(uMinor, vMinor)];
//             F.addEdge(uG0, vG0, ew);
//         });
//     };

//     auto retraction = [&](bool force = false){
//         assert(!force || Cz.size() > 1);
//         vector<edge> outsideEdges;
//         Graph GMinor;
//         map<NodePair, NodePair> edge_GMinor_to_G0;
//         map<node, node> node_GT_to_GMinor;
//         auto& CzLast = Cz[Cz.size() - 1];

//         for (node uGT: CzLast) {
//             node uGMinor = GMinor.addNode();
//             node_GT_to_GMinor[uGT] = uGMinor;
//         } 

//         for (node uGT: CzLast) {
//             GT.forEdgesOf(uGT, [&](node uGT2, node vGT, edgeweight ew) {
//                 assert(uGT == uGT2);
//                 if (CzLast.contains(vGT)) {
//                     auto uGMinor = node_GT_to_GMinor[uGT];
//                     auto vGMinor = node_GT_to_GMinor[vGT];
//                     GMinor.addEdge(uGMinor, vGMinor, ew);

//                     auto uvG0 = edge_GT_to_G0[minmax(uGT, vGT)];
//                     edge_GMinor_to_G0[minmax(uGMinor, vGMinor)] = uvG0;
//                 } else {
//                     outsideEdges.push_back(edge{
//                         uGT,
//                         vGT,
//                         ew,
//                         ew
//                     });
//                 }
//             });
//         }
//         minorMst(GMinor, edge_GMinor_to_G0);

//         //
//         // Now need to change the GT graph
//         //  This will also update GT maps...
//         //
//         for (node uGT: CzLast) {
//             removeGTNode(uGT);
//         }
//         node contractedNode = GT.addNode();

//         map<node, edgeweight> smallestEdge;
//         for (edge e: outsideEdges) {
//             assert(CzLast.contains(e.u) && !CzLast.contains(e.v));
//             if (!smallestEdge.contains(e.v) || e.key < smallestEdge[e.v]) {
//                 smallestEdge[e.v] = e.key;
//             }
//         }

//         for(auto [vGT, ew]: smallestEdge) {
//             GT.addEdge(contractedNode, vGT, ew);
//         }

//         updateMinLinks(minLink.size() - 1);
//         chainLink.pop_back();
//         Cz.pop_back();
//         // [TODO heaps]
//         // updateHeaps(heaps.size() - 1);
//     };

//     auto needFusion = [&](){

//         return false;
//     };
//     auto heapsEmpty = [&](){
//         // for (auto& hs : heaps) {
//         //     for (auto& h : hs) {
//         //         if (!h.empty()) {
//         //             return false;
//         //         }
//         //     }
//         // }
//         // return true;
//         return heaps.empty();
//     };
//     auto shouldRetract = [&](){return (rand() & 0xFF) == 0  && Cz.size() > 1;};

//     while (!heapsEmpty()) {
//         if (shouldRetract()) {
//             retraction();
//         } else {
//             edge e = fusion();
//             // if (needFusion()){
//             //     fusion();
//             // }
//             extension(e);
//         }
//     }
//     while (!Cz.empty()) {
//         retraction(true);
//     }
//     // [TODO] do retractions to all leftover T hierarchy
//     // while (Cz.size() > 0) { retraction(true); } // true because we want to remove last Cz... 

//     // [TODO] REMOVE THIS
//     G0.forEdges([&](node u, node v, edgeweight ew){
//         if(rand() & 1) {
//             F.addEdge(u, v, ew);
//         } else {
//             B.addEdge(u, v, ew);
//         }
//     });

//     // [STEP 4]
//     // This step is done during [STEP 3] retractions.

//     // [STEP 5]
//     B.forEdges([&](node u, node v, edgeweight ew) {
//         F.addEdge(u, v, ew);
//     });
//     Graph res = mst(F, t);
//     res.forEdges([&](node u, node v, edgeweight ew) {
//         auto [uG, vG] = edge_G0_to_G[minmax(u, v)];
//         forest.addEdge(uG, vG, ew);
//     });
//     return forest;
// }

// Max number of nodes ~ billion
std::pair<int, std::vector<int>> tHierarchySize(int n, int m) {
    if (n < 3) {
        return {{3}, {10,4,1}};
    }
    // Arbitrary sequence for now that grows really fast -- no point in using Ackermann sequence...
    // int desiredNumberOfLeaves[7] = {1, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 30, 1 << 30};
    int desiredNumberOfLeaves[7] = {1, 1 << 2, 1 << 6, 1 << 16, 1 << 30, 1 << 30, 1 << 30};
    // int desiredNumberOfLeaves[7] = {1,  1 << 2, 1 << 30, 1 << 30, 1 << 30, 1 << 30, 1 << 30};
    int desiredNumberOfChildren[7];
    desiredNumberOfChildren[0] = 1;
    for (int i = 1; i < 7; ++i) desiredNumberOfChildren[i] = desiredNumberOfLeaves[i] / desiredNumberOfLeaves[i - 1];

    std::vector<int> desiredSize;
    int d = 0;
    for (; desiredNumberOfLeaves[d] < n; ++d) {}
    for (int i = d; i >= 0; --i) desiredSize.push_back(desiredNumberOfChildren[i]);

    return {d, desiredSize};
}

NetworKit::Graph Chazelle2000MinimumSpanningTree::msf(NetworKit::Graph G, int t) {
    // std::cout << "Inside msf" << std::endl;
    auto tree = NetworKit::GraphTools::copyNodes(G);
    auto connected_components = NetworKit::ConnectedComponents(G);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (auto i = 0; i < connected_components.numberOfComponents(); i++) {
        // std::cout << "Inside cc " << i << std::endl;
        auto subgraph = NetworKit::GraphTools::subgraphFromNodes(G, [](const auto& components){
            std::unordered_set<NetworKit::node> s;
            for (auto n: components){
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
            // std::cout << u << " <-> " << n << std::endl; 
        });
        subgraph = NetworKit::GraphTools::getCompactedGraph(subgraph, mappingTo);
        // subgraph.forNodes([](node u){std::cout<< u << " ";}); std::cout << std::endl;
        // subgraph.forEdges([](node u, node v, edgeweight ew){std::cout << "{"<< u << " " << v << " "<< ew << "} ";}); std::cout << std::endl;

        // std::cout << "Calling mst on " << i << std::endl;
        auto tree1 = mst(subgraph, t);
        tree1.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
            // std::cout << "{"<< u << " " << v << " "<< ew << "} ";
            tree.addEdge(mappingFrom[u], mappingFrom[v], ew);
            // std::cout << "["<< mappingFrom[u] << " " << mappingFrom[v] << " "<< ew << "] ";
        });
        // std::cout << std::endl;
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
        s.insert(std::minmax(u,v));
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

    // std::cout << "////////////////////////////////\n";
    // G.forEdges([](node u, node v, edgeweight ew){
    //     std::cout << u << " " << v << " " << ew << '\n';
    // });
    // G.forNodes([](node u){
    //     std::cout << u << " ";
    // });
    // std::cout << '\n';
    // std::cout << "////////////////////////////////" << std::endl;
    
    // std::cout << "INSIDE CHAZELLE" << std::endl;
    // [STEP 1]
    if (t <= 1 || G.numberOfNodes() < MIN_NUMBER_NODES) {
        BoruvkaMinimumSpanningTree boruvka(G);
        boruvka.run();
        return boruvka.getForest();
    }
    
    // [STEP 2]
    std::cout << "BEFORE " << C << " BORUVKA STEPS" << std::endl;
    auto [G0, edge_G0_to_G, forest] = boruvkaSteps(G, C);
    std::cout << "G0 SIZE " << G0.numberOfNodes() << " " << G0.numberOfEdges() << std::endl;
    if (G0.numberOfEdges() == 0) {
        std::cout << "BORUVKA STEPS FINISHED" << std::endl;
        // assert(false);
        return forest;
    }

    std::cout << "////////////////////////////////\n";
    // G0.forEdges([](node u, node v, edgeweight ew){
    //     std::cout << u << " " << v << " " << ew << '\n';
    // });
    // G0.forNodes([](node u){
    //     std::cout << u << " ";
    // });
    std::cout << '\n';
    std::cout << "////////////////////////////////" << std::endl;

    // edge dummy_edge{0, 0, 0, 0, 0, false, false};
    // SoftHeap<edge*> heaps(0.1);
    vector<vector<SoftHeap<edge*>>> heaps;
    // heaps.extractMin();
    // std::priority_queue<edge> heaps;
    set<NodePair> contractedEdges;
    set<NodePair> badEdges;
    vector<edge> edges(G0.numberOfEdges());
    // std::list<edge*> justCorruptedElementsList;
    vector<set<int>> Cz;
    vector<int> parent(G0.numberOfNodes());
    vector<bool> fusionNode(G0.numberOfNodes(), false);
    vector<vector<pair<NodePair, double>>> minLink;
    map<NodePair, int> G0_edge_id;
    vector<bool> visited(G0.numberOfNodes(), false);
    assert(!G.isDirected());
    assert(!G0.isDirected());
    auto [d, desiredSize] = tHierarchySize(G0.numberOfNodes(), G0.numberOfEdges());


    // we have Cz0, Cz1, .... Czk
    auto k = [&](){return Cz.size() - 1;};

    // [TODO] Check all edge{} constructions AFTER THIS and change to this...
    auto makeEdge = [&](node u, node v, edgeweight ew) {
        return edge{u, v, ew, ew, G0_edge_id[minmax(u, v)], false, false};
    };

    auto minBorderEdge = [&](){
        double minKey = __DBL_MAX__;
        int mini = -1, minj = -1;

        // edge minEdge = heaps.top();
        // heaps.pop();
        // while(edges[minEdge.id].removed) {
        //     minEdge = heaps.top();
        //     heaps.pop();
        // }
        // return std::tuple{minEdge, 0, 0};
        // [TODO] Implement heaps
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
        // while(minEdge->removed) {
            // minEdge = heaps.extractMin();
            // if (minEdge->corrupted) {
            //     badEdges.insert(minmax(minEdge->u, minEdge->v));
            // }
        // }
        
        return std::tuple{*minEdge, 0, 0};
    };

    auto findAndDeleteEdgeFromHeaps = [&](node u, node v) {
        // [TODO] Implement heaps
        auto eid = G0_edge_id[minmax(u, v)];
        edges[eid].removed = true;
        if (edges[eid].corrupted) {
            // [TODO] When heaps are implemented it can be moved to retraction and fusion... - inside the updateHeaps or something like that
            badEdges.insert(minmax(u, v));
        }
        return edges[eid];
    };

    auto insertNewBorderEdge = [&](edge* e) {
        heaps[k()][k()] = std::move(insert(std::move(heaps[k()][k()]), e));
    };

    auto insertNewBorderEdges = [&](const vector<edge>& newBorderEdges) {
        // [TODO] Implement heaps
        for (auto e: newBorderEdges) {
            NodePair np = minmax(e.u, e.v);
            auto eid = G0_edge_id[np];
            edges[eid] = e;
            // [TODO] Implement heaps
            insertNewBorderEdge(&edges[eid]);
            // heaps.push(e);
        }
        return;
    };
    
    auto nodeTopParent = [&](node u) {
        int p = static_cast<int>(u);
        while(parent[p] != p) p = parent[p];
        return p;
    };

    auto leftmostSmallerMinLink = [&](int ckey){
        // int i = 0;
        // int j = 0;
        NodePair ab;

        // std::cout << "*************\n";
        // for (int i = 0; i < k(); ++i) {
        //     for (int j = i + 1; j <= k(); ++j) {
        //         std::cout << minLink[i][j].second << " ";
        //     }std::cout << '\n';
        // }
        for (int i = 0; i < k(); ++i) {
            for (int j = i + 1; j <= k(); ++j) {
                // std::cout << minLink[i][j].second << " ";
                // doesnt exist
                if (minLink[i][j].second < 0) continue;
                if (minLink[i][j].second <= ckey) {
                    return tuple<int, int, NodePair>{i, j, minLink[i][j].first};
                }
            }
            // std::cout << '\n';
        }
        std::cout << "*************\n";
        std::cout << "*************\n";
        
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
                // if (minLink[i][j].first.first != -1) {
                //     assert(edges[G0_edge_id[minLink[i][j].first]].ckey == minLink[i][j].second);
                // }
                if (i >= k() || j > k()) {
                    minLink[i][j] = {{-1, -1}, __DBL_MAX__};
                }
            }
        }

        // for (auto e: insideEdges) {
        //     int u = e.u, v = e.v;
        //     int up = nodeTopParent(u), vp = nodeTopParent(v);
        //     int ui = -1, vi = -1;
        //     for (int i = 0; i <= k(); ++i) {
        //         if (Cz[i].contains(up)) ui = i;
        //         if (Cz[i].contains(vp)) vi = i;
        //     }
        //     assert(ui != -1);
        //     assert(vi != -1);
        //     if (ui > vi) {
        //         swap(ui, vi);
        //     }
        //     if (ui != vi) {
        //         assert(minLink[ui][vi].second <= e.key);
        //         assert(minLink[ui][vi].second >= 0);
        //     }
        // }
    };

    auto updateHeaps = [&](int last_i) {
        // [TODO] when heaps
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

    auto popHeaps = [&]() {
        heaps[k() - 1][k() - 1] = std::move(meld(std::move(heaps[k()][k()]), std::move(heaps[k() - 1][k() - 1])));
        heaps.pop_back();
        for (auto& h: heaps) {
            h.pop_back();
        }
    };

    auto retraction = [&]() {
        assert(Cz.size() > 1);
        int newP = parent.size();
        parent.push_back(newP);
        fusionNode.push_back(false);
        auto vertices = Cz[Cz.size() - 1];
        popHeaps();
        // heaps[k() - 1][k() - 1] = std::move(meld(std::move(heaps[k()]), std::move(heaps[k() - 1])));
        // heaps.pop_back();
        Cz.pop_back();
        if (vertices.size() == 1) {
            parent.pop_back();
            fusionNode.pop_back();
            newP = *vertices.begin();
        } else {
            for (auto v: vertices) {
                parent[v] = newP;
            }
        }
        
        Cz[Cz.size() - 1].insert(newP);
        
        updateMinLinks();
    };

    auto fusion = [&]() -> edge {
        auto [uvEdge, mini, minj] = minBorderEdge();
        node uG0{uvEdge.u}, vG0{uvEdge.v};
        auto [link_i, link_j, ab] = leftmostSmallerMinLink(uvEdge.ckey); // [TODO CHECK] ckey ???
        // std::cout << "Got the edge!" << std::endl;
        // std::cout<<link_i<<" "<<link_j<<" "<<ab.first<<" "<<ab.second<<std::endl;
        // std::cout<<uvEdge.u<<" "<<uvEdge.v << " " <<uvEdge.key <<std::endl;
        if (link_i == -1) return uvEdge;
        // while (link_i < k()) {
        //     std::cout<<k()<<'\n';
        //     retraction();
        // }
        // return uvEdge;
        
        auto [a, b] = ab;
        auto ap = nodeTopParent(a); 
        auto bp = nodeTopParent(b);
        if (!Cz[link_i].contains(ap)) {
            swap(a, b);
            swap(ap, bp);
        }
        assert(Cz[link_i].contains(ap));

        int newP = parent.size();
        parent.push_back(newP);
        fusionNode.push_back(true);
        while (k() > link_i) {
            for (auto v: Cz[k()]) {
                parent[v] = newP;
            }
            popHeaps();
            // heaps[k() - 1] = std::move(meld(std::move(heaps[k()]), std::move(heaps[k() - 1])));
            // heaps.pop_back();
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
        // if (minLink[i][j].second > e.key || minLink[i][j].second < 0) {
        //     minLink[i][j] = {minmax(e.u, e.v), e.key};
        // }
        // if (minLink[i][j].second < 0) minLink[i][j].first = {-1, -1};
    };

    int extensions = 0;
    auto extension = [&](edge e) {
        // std::cout << "BRUH 2" <<std::endl;
        extensions += 1;
        node u = e.u, v = e.v;
        if (!visited[u]) swap(u, v);
        assert(visited[u]);
        assert(!visited[v]);
        // std::cout << "BRUH 3" <<std::endl;

        heaps.push_back(vector<SoftHeap<edge*>>());
        for (int i = 0; i < Cz.size(); ++i) heaps[Cz.size()].push_back(SoftHeap<edge*>(0.1));
        Cz.push_back({v});
        for (int i = 0; i < Cz.size(); ++i) heaps[i].push_back(SoftHeap<edge*>(0.1));
        visited[v] = true;
        // std::cout << "BRUH 4" <<std::endl;
        G0.forEdgesOf(v, [&](node vv, node w, edgeweight ew) {
            // std::cout << "BRUH 5" <<std::endl;
            assert(v == vv);
            if (visited[w]) {
                edge e = findAndDeleteEdgeFromHeaps(v, w);
                updateMinLinksWithEdge(e);
                insideEdges.insert(e);
            } else {
                int eid = G0_edge_id[minmax({v, w})];
                // [TODO] impelment heaps
                heaps[k()][k()] = std::move(insert(std::move(heaps[k()][k()]), &edges[eid]));
                // heaps.push(edges[eid]);
            }
        });
        // std::cout << "BRUH 6" <<std::endl;
        updateMinLinks();
        // std::cout << "BRUH 7" <<std::endl;
    };

    auto shouldRetract = [&]() {
        return Cz.size() == d || (Cz[k()].size() >= desiredSize[k()] && k() > 0);
        // return (rand() & 0xFF) == 0  && Cz.size() > 1;
    };

    auto shouldFinish = [&]() {
        return extensions >= G0.numberOfNodes() - 1;
    };

    auto initialization = [&]() {
        std::cout << "Inside initialization" << std::endl;
        int eidNow = 0;
        G0.forEdges([&](node u, node v, edgeweight ew){
            edges[eidNow] = edge{u, v, ew, ew, eidNow, false, false};
            G0_edge_id[minmax(u, v)] = eidNow;
            eidNow += 1;
            // std::cout << "First edge loop" << std::endl;
        });

        // G0.forNodes([](node u){
        //     std::cout << u << '\n';
        // });
        
        // std:: cout << "0 node degree: " << G0.degree(0) << '\n';
        // std::cout << "Inside initialization 2" << std::endl;
        for(int i = 0; i < parent.size(); i++) parent[i] = i;
        Cz.push_back({0});
        heaps.push_back({SoftHeap<edge*>(0.1)});
        G0.forEdgesOf(0, [&](node u, node v, edgeweight ew) {
            // std::cout << "Second edge loop" << std::endl;
            int eid = G0_edge_id[minmax(u, v)];
            // std::cout << u << ' ' << v << ' ' << ew << ' ' << eid << ' ' << edges[eid].toString() << std::endl;
            // [TODO] implement heaps
            heaps[0][0] = std::move(insert(std::move(heaps[0][0]), &edges[eid]));
            // heaps.push(edges[eid]);
        });
        std::cout << "Inside initialization 3" << std::endl;
        minLink = vector<vector<pair<NodePair, double>>>(d, vector<pair<NodePair, double>>(d, {{-1, -1}, __DBL_MAX__}));
        visited[0] = true;
        // SoftHeap<edge*>::justCorruptedElementsList = &justCorruptedElementsList;
    };

    std::cout << "Before initialization" << std::endl;
    initialization();
    std::cout << "After initialization" << std::endl;
    while (!shouldFinish()) {
        // for (auto cz : Cz) {
        //     for (auto v : cz) {
        //         std :: cout << v << " "; 
        //     }std::cout << '\n';
        // }std::cout << "--------\n";
        if (shouldRetract()) {
            // std:: cout << "retraction " << G0.numberOfNodes() - 1 - extensions << std::endl; 
            retraction();
        } else {
            // std:: cout << "fusion + extension " << G0.numberOfNodes() - 1 - extensions << std::endl;
            edge e = fusion();
            // if (needFusion()){
            //     fusion();
            // }
            extension(e);
        }
    }
    std::cout << "Processed the whole graph!" << std::endl;
    // [TODO] consider Cz.size() == 1 carefully
    while (Cz.size() > 1) {
        int newP = parent.size();
        parent.push_back(newP);
        fusionNode.push_back(false);
        auto vertices = Cz[Cz.size() - 1];

        heaps[k() - 1][k() - 1] = std::move(meld(std::move(heaps[k()][k()]), std::move(heaps[k() - 1][k() - 1])));
        heaps.pop_back();
        Cz.pop_back();

        for (auto v: vertices) {
            parent[v] = newP;
        }
        Cz[Cz.size() - 1].insert(newP);
    }
    int root = parent.size();
    parent.push_back(root);
    fusionNode.push_back(false);
    for (int v: Cz[0]) {
        parent[v] = root;
    }
    assert(parent.size() == fusionNode.size());

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
        assert(v == root);

        int d = depth[i];
        v = i;
        while (v != parent[v]) {
            int prev_v = v;
            v = parent[v];
            d -= 1;

            depth[v] = d;
            int new_v = fromG0Maps[v].size();
            // std::cout << v << ' ' << fromG0Maps[v].size() << std::endl;
            fromG0Maps[v].insert({prev_v, new_v});
        }
    }
    std :: cout << "depth calculated" << std::endl;
    
    vector<Graph> CzGraphs(parent.size(), Graph(0, true));
    for (int i = G0.numberOfNodes(); i < parent.size(); ++i) {
        CzGraphs[i] = Graph(fromG0Maps[i].size(), true);
        // for (int j = 0; j < fromG0Maps[i].size(); ++j) CzGraphs[i].addNode();
    }
    std :: cout << "Cz graphs calculated" << std::endl;

    
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

        while(l != r) {
            prev_l = l; prev_r = r;
            l = parent[l]; r = parent[r];
        }

        return tuple<int, int, int>{l, prev_l, prev_r};
    };

    vector<map<NodePair, double>> smallestEdge(CzGraphs.size());

    G0.forEdges([&](node u, node v, edgeweight ew) {
        if (badEdges.contains(minmax(u, v))) return;

        
        auto [lca, l, r] = lcaLR(u, v);
        // for (auto [key, val]: fromG0Maps[lca]) {
        //     std :: cout << "{ "<< key << " -> " << val << "}, ";
        // }std :: cout << std::endl;
        auto minorL = fromG0Maps[lca][l];
        auto minorR = fromG0Maps[lca][r];
        // std :: cout << "lca " << lca << " l " << l << " r " << r << " minorL " << minorL << " minorR " << minorR << std::endl; 
        // std :: cout << "N: " << CzGraphs[lca].numberOfNodes() << " " << fromG0Maps[lca].size() << std::endl;
        // CzGraphs[lca].addEdge(minorL, minorR, ew);
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
        // std :: cout << "adding edges " << i << std::endl;
        for (auto [key, val]: smallestEdge[i]) {
            CzGraphs[i].addEdge(key.first, key.second, val);
        }
    }

    std :: cout << "Cz graphs edges added" << std::endl;

    Graph G0copy = G0;
    BoruvkaMinimumSpanningTree bbb(G0copy);
    bbb.run();
    auto actualMst = bbb.getForest();
    set<NodePair> actualMstEdges;
    actualMst.forEdges([&](node u, node v, edgeweight ew) {
        actualMstEdges.insert(minmax(u, v));
    });

    // auto assertContractible = [&](int i) {
    //     std::cout <<"Assert " << i << " contractible\n";

    //     CzGraphs[i].forEdges([&](node u, node v, edgeweight ew){
    //         std :: cout << "Assert" << actualMstEdges.contains(toG0edgeMaps[i][minmax(u, v)]) << '\n';
    //         assert(actualMstEdges.contains(toG0edgeMaps[i][minmax(u, v)]));
    //     });
    //     std::cout << "GREPME\n";
    // };

    std::cout << " B EDGES SIZE " << badEdges.size() << '\n';
    // assert(badEdges.size() == 0);

    for (int i = G0.numberOfNodes(); i < CzGraphs.size(); ++i) {
        // std::cout << "recursive call depth: " << depth[i] << " size: " << CzGraphs[i].numberOfNodes() << std::endl;
        // std::cout << i << " out of " << CzGraphs.size() << std::endl;
        // CzGraphs[i].forNodes([](node u){std::cout<< u << " ";}); std::cout << std::endl;
        // CzGraphs[i].forEdges([](node u, node v, edgeweight ew){std::cout << "{"<< u << " " << v << " "<< ew << "} ";}); std::cout << std::endl;
        if (fusionNode[i]) continue;
        CzGraphs[i] = msf(CzGraphs[i], t - 1);
        // assertContractible(i);
    }
    std :: cout << "msts called" << std::endl;
    
    auto F = NetworKit::GraphTools::copyNodes(G0);
    for (int i = G0.numberOfNodes(); i < CzGraphs.size(); ++i) {
        // std ::cout << "Processing CzGraphs[" << i << "]" << std::endl;
        CzGraphs[i].forEdges([&](node u, node v){
            auto [x, y] = toG0edgeMaps[i][minmax(u, v)];
            F.addEdge(x, y, edges[G0_edge_id[minmax(x, y)]].key);
            // std::cout << x << " " << y << std::endl; 
        });
    }
    std :: cout << "edges from step 4 added" << std::endl;
    for (auto [x, y]: badEdges) {
        F.addEdge(x, y, edges[G0_edge_id[minmax(x, y)]].key);
    }
    std :: cout << "bad edges added" << std::endl;

    // [STEP 5]
    Graph res = mst(F, t);
    // forest.forEdges([&](node u, node v, edgeweight ew) {
    //     std :: cout << "forest edge: " << u << " " <<  v << " " << ew << '\n';
    // });

    res.forEdges([&](node u, node v, edgeweight ew) {
        // std :: cout << "res edge G0: " << u << " " <<  v << " " << ew << std::endl;
        auto [uG, vG] = edge_G0_to_G[minmax(u, v)];
        forest.addEdge(uG, vG, ew);
        // std :: cout << "res edge  G: " << uG << " " <<  vG << " " << ew << std::endl;
    });

    for(int i = 0; i < parent.size(); ++i) {
        int v = i;
        // std :: cout << v;
        while(v != parent[v]) {
            // std :: cout << " -> " << parent[v];
            v = parent[v]; 
        }
        // std ::cout << std::endl;
    }

    return forest;
}

}  /* namespace Koala */
