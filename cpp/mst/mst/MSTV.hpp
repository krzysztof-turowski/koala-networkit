#pragma once

#include "MinimalSpanningTree.hpp"
#include "BoruvkaMST.hpp"

namespace MST {
    /*
     * O(V+E) Minimal Spanning Tree Verification Algorithm. Returns true if a given mst is truly minimal.
     * Torben Hagerup, 2010: "An Even Simpler Linear-Time Algorithm for Verifying Minimum Spanning Trees"
     */
    bool HagerupMSTV(const MinimalSpanningTree& mst);
    bool HagerupMSTV(const Graph &mst, const Graph &G);

    /*
     * Similar to HagerupMSTV.
     * F is a spanning tree of G.
     * For each edge e=(u,v) in G and a path P=<u, ..., v> in F, this function removes e if w(e) > w(e') for e' in P.
     */
    void removeF_HeavyEdges(Graph& G, const Graph& F);
    void removeF_HeavyEdges(Graph& G, const Graph& F, Graph& G2, const Graph& F2, const std::vector<node> &S_invert);
}

#include <structures/LCA.hpp>
#include "MinimalSpanningTree.hpp"
#include <iostream>

void checkMediansCorrectness(const std::vector<MST::count>& medians) {
    for(int i = 0; i<medians.size(); i++) {
        if (i == 0) {
            assert(medians[i] == 0);
            continue;
        }
        std::vector<int> elements;
        auto hset = i;
        for(int j = 0; (1 << j) <= hset; j++) {
            if ((1<<j) & hset) {
                elements.push_back(j);
            }
        }
        assert(elements.size() > 0);
        assert(elements.at(elements.size()/2) == medians[i]);
    }
    // std::cout <<"Medians correct!" << std::endl;
}

namespace MST {
    using h_set = uint64_t; // bitset of depths. ith bit corresponds to ith depth.

    template <typename FBT>
    edgeweight edgeWeightToParent(node u, const FBT& fbt) {
        // value for root not explained in the paper.
        if (u == fbt.getRoot()) return 0;
        auto w = fbt.getTree().weight(fbt.getParent(u).value(), u);
        assert(w > 0); // found appropriate edge
        return w;
    };

    template <typename FBT>
    std::vector<node> treePathMaxima(
        FBT& fbt, const std::vector<node>& upper, const std::vector<node>& lower) {

        /// This code is derived from the paper which has an implementation in D.
        /// It resembles the original D code as much as possible - including multiple lambdas instead of
        /// top-level functions.
        /// There was a bug in the paper due to operators order - see `visit`!

        count height = 0; // height of `fbt` = maximum depth
        count n = fbt.getTree().numberOfNodes(), m = upper.size();
        // depth: node -> depth.
        std::vector<count> depth(n);
        // median: h_set -> element of h_set. Example: 0b1011 represents {0,1,3}. median[0b1011] = 1;
        std::vector<index> median;
        // `L` and `Lnext` combined form a linkedlist of "indices of queries" (query_ids).
        // See the code where it's populated for more explanation.
        std::vector<index> L(n, NetworKit::none), Lnext(m, NetworKit::none);
        // answer: query_id -> node. This vector is returned from the function.
        std::vector<node> answer(m, NetworKit::none);
        // P: depth -> node. DFS stack, used in `visit`.
        std::vector<node> P;
        // D[u]: a set of depths of endpoints above u of query paths that contain u. D[root] = 0.
        std::vector<h_set> D(n, 0);

        auto init = [&](node u, count u_depth, auto&& init) -> void {
            depth[u] = u_depth;
            if (u_depth > height)
                height = u_depth;
            for (auto i = L[u]; i != NetworKit::none; i = Lnext[i])
                D[u] |= 1 << depth[upper[i]]; // this is executed only for leaves of `fbt`.

            fbt.getTree().forNeighborsOf(u, [&](NetworKit::node, NetworKit::node child, NetworKit::edgeweight) {
                init(child, u_depth + 1, init);
                // exclude `u`. no-op for leaves and works recursively up to the root.
                D[u] |= (D[child] & ~(1 << u_depth));
            });
        };

        auto median_table = [&median](count h) -> void {
            // Fills a table of size 2^(h+1) whose entry in position i, for
            // i=0,...,2^(h-1)-1, is the median of the set represented by i.

            std::vector<h_set> T((1 << h) + 1);
            median.resize(1 << (h + 1));

            auto subsets = [&](count n, count k, index p, auto&& subsets) -> index {
                // Stores the subsets of size k of {0,...,n-1} in T,

                // starting in position p, and returns p plus their number.
                if (n < k)
                    return p;
                if (k == 0) {
                    T[p] = 0;
                    return p + 1;
                }

                index q = subsets(n - 1, k - 1, p, subsets);
                for (auto i = p; i < q; i++)
                    T[i] |= 1 << (n - 1);
                return subsets(n - 1, k, q, subsets);
            };

            for (count s = 0; s <= h; s++)
                for (count k = 0; k <= s; k++) {

                    auto p = subsets(h - s, k, 0, subsets);
                    auto q = subsets(s, k, p, subsets);
                    q = subsets(s, k + 1, q, subsets);
                    for (count i = 0; i < p; i++) {
                        auto b = (1 << (s + 1)) * T[i] + (1 << s); // fixed high bits
                        for (auto j = p; j < q; j++)
                            median[b + T[j]] = s; // variable low bits
                    }
                }
            checkMediansCorrectness(median);
            return;
        }; // end median_table

        auto down = [](h_set A, h_set B) -> h_set {
            // Returns A "downarrow" B
            return B & (~(A | B) ^ (A + (A | ~B)));
        };

        auto visit = [&](node v, h_set S, auto&& visit) -> void {
            // when called, S is "S of `parent(v)`" or \emptyset for root.

            auto binary_search = [&](double w, int S) -> count {
                // Returns max({j in S | weight[P[j]]>w} union {0})
                // needed for Sv definition on the bottom of paper's 183 page.

                if (S == 0)
                    return 0;
                auto j = median[S];
                for(; S != (1 << j); j = median[S]) { // `while |S|>1` (or, to be more specific, `while S != {j}`).
                    S &= (edgeWeightToParent(P[j], fbt) > w) ? (~((1 << j) - 1)) : ((1 << j) - 1);
                }
                return (edgeWeightToParent(P[j], fbt) > w) ? j : 0; // if not found, return `-inf` (or 0).
            };


            P[depth[v]] = v; // push current node on stack
            // sup{j' \in down(Dv, Su) : w(Pv(j')) > w(v)}
            int k = binary_search(edgeWeightToParent(v, fbt), down(D[v], S));

            S = down(D[v], S & ((1 << (k + 1)) - 1) | (1 << depth[v]));
//            S = down(D[v], S & (1 << (k + 1) - 1) | (1 << depth[v])); - it was a BUG in the paper!
            for (auto i = L[v]; i != NetworKit::none; i = Lnext[i]) {
                auto singleton_depth_of_ans = down(1 << depth[upper[i]], S);
                auto depth_of_ans = median[singleton_depth_of_ans];
                answer[i] = P[depth_of_ans];
            }

            fbt.getTree().forNeighborsOf(v, [&](NetworKit::node, NetworKit::node child, NetworKit::edgeweight) {
                visit(child, S, visit);
            });
        }; // end visit

        for (auto i = 0; i < m; i++) { // distribute queries to lower nodes.
            // L[u] - beginning of a linked list query_ids that have lower node `u`.
            // L[u] ==  index of the first query. Lnext[L[u]] - index of the next query. And so on until we hit `Networkit::none`.
            // L[] is indexed by nodes and Lnext[] is indexed by query_ids.
            // L[u] points to the beginning of the linked list for `u`. This list is is stored in `Lnext[]`.
            Lnext[i] = L[lower[i]];
            L[lower[i]] = i;
        }

        init(fbt.getRoot(), 0, init);
        P.resize(height + 1);
        median_table(height);
        visit(fbt.getRoot(), 0, visit);
        return answer;
    } // end tree_path_maxima

    bool HagerupMSTV(const Graph &mst, const Graph &G) {
        std::cout << "CHECK" << std::endl;
        // Boruvka runs in linear time on trees.
        BoruvkaMST fbtOnMst(mst, true, true);

        std::vector<std::tuple<NetworKit::node, NetworKit::node, NetworKit::edgeweight>> G_minus_M;
        G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
             if (!mst.hasEdge(u, v)) {
                G_minus_M.push_back({u, v, w});
             }
        });
        Koala::LCA<FullBranchingTree> lca(fbtOnMst.fullBranchingTree.value());
        std::vector<node> upper, lower;
        for (auto [u, v, w] : G_minus_M) {
            auto lowestCommonAncestor = lca.query(u, v);
            upper.push_back(lowestCommonAncestor), upper.push_back(lowestCommonAncestor);
            lower.push_back(u), lower.push_back(v);
        }

        auto answers = treePathMaxima(fbtOnMst.fullBranchingTree.value(), upper, lower);
        for (index i = 0; i < answers.size(); i++) {
            if (edgeWeightToParent(answers[i], fbtOnMst.fullBranchingTree.value()) > std::get<2>(G_minus_M[i / 2])) {
                return false;
            }
        }
        return true;
    }

    bool HagerupMSTV(const MinimalSpanningTree& mst, Graph &G) {
        if (!mst.getOriginalGraph().has_value()) {
            throw std::invalid_argument("Minimal Spanning Tree verification needs access to the original graph!");
        }
        // Boruvka runs in linear time on trees.
        BoruvkaMST fbtOnMst(mst.getSpanningTree(), true, true);

        auto G_minus_M = mst.getEdgesNotBelongingToMST();
        Koala::LCA<FullBranchingTree> lca(fbtOnMst.fullBranchingTree.value());
        std::vector<node> upper, lower;
        for (auto e: G_minus_M) {
            node lowestCommonAncestor = lca.query(e.first.first, e.first.second);
            upper.push_back(lowestCommonAncestor);
            lower.push_back(e.first.first);
            upper.push_back(lowestCommonAncestor);
            lower.push_back(e.first.second);
        }

        auto answers = treePathMaxima(fbtOnMst.fullBranchingTree.value(), upper, lower);
        for (index i = 0; i < answers.size(); i++) {
            if (edgeWeightToParent(answers[i], fbtOnMst.fullBranchingTree.value()) > G_minus_M[i / 2].second) {
                return false;
            }
        }
        return true;
    }

    void removeF_HeavyEdges(Graph& G, const Graph& F) {
        /// This algorithm is very similar to HagerupMSTV. However, HagerupMSTV terminates whenever it finds an
        /// F_heavy edge - hence we cannot reuse it.
        BoruvkaMST FbtOnF(F, true, true);
        Koala::LCA<FullBranchingTree> lca(FbtOnF.fullBranchingTree.value());
        std::vector<node> upper, lower;
        std::vector<std::pair<NodePair, edgeweight>>  edges;
        G.forEdges([&](const node u, const node v, edgeweight ew) {
           node lowestCommonAncestor = lca.query(u, v);
           upper.push_back(lowestCommonAncestor);
           lower.push_back(u);
           upper.push_back(lowestCommonAncestor);
           lower.push_back(v);
           edges.push_back({{u, v}, ew});
        });
        auto answers = treePathMaxima(FbtOnF.fullBranchingTree.value(), upper, lower);
        for (index i = 0; i < answers.size(); i+=2) {
            node u = edges[i/2].first.first;
            node v = edges[i/2].first.second;
            edgeweight w = edges[i / 2].second;
            edgeweight maxOnThe1stHalf = edgeWeightToParent(answers[i], FbtOnF.fullBranchingTree.value());
            edgeweight maxOnThe2ndHalf = edgeWeightToParent(answers[i+1], FbtOnF.fullBranchingTree.value());
            if (w > maxOnThe1stHalf && w > maxOnThe2ndHalf) {
                G.removeEdge(u, v);
            }
        }
    }

    void removeF_HeavyEdges(Graph& G, const Graph& F, Graph& G2, const Graph& F2, const std::vector<node> &S_invert) {
        /// This algorithm is very similar to HagerupMSTV. However, HagerupMSTV terminates whenever it finds an
        /// F_heavy edge - hence we cannot reuse it.
        BoruvkaMST FbtOnF(F, true, true);
        Koala::LCA<FullBranchingTree> lca(FbtOnF.fullBranchingTree.value());
        std::vector<node> upper, lower;
        std::vector<std::pair<NodePair, edgeweight>>  edges;
        G.forEdges([&](const node u, const node v, edgeweight ew) {
           node lowestCommonAncestor = lca.query(u, v);
           upper.push_back(lowestCommonAncestor), upper.push_back(lowestCommonAncestor);
           lower.push_back(u), lower.push_back(v);
           edges.push_back({{u, v}, ew});
        });
        auto answers = treePathMaxima(FbtOnF.fullBranchingTree.value(), upper, lower);
        for (index i = 0; i < answers.size(); i+=2) {
            auto [u, v] = edges[i/2].first;
            edgeweight w = edges[i / 2].second;
            edgeweight maxOnThe1stHalf = edgeWeightToParent(answers[i], FbtOnF.fullBranchingTree.value());
            edgeweight maxOnThe2ndHalf = edgeWeightToParent(answers[i+1], FbtOnF.fullBranchingTree.value());
            if (w > maxOnThe1stHalf && w > maxOnThe2ndHalf) {
                G.removeEdge(u, v);
                G2.removeEdge(S_invert[u], S_invert[v]);
            }
        }
    }
}
