#include <bits/stdc++.h>

#include <boost/functional/hash.hpp>

#include <networkit/graph/GraphTools.hpp>

#include "krt_edge_designator.hpp"

#include "krtgraph.hpp"
#include "logger.hpp"

#include "dynamic_trees/dyn_tree.h"

#ifndef PII
#define PII std::pair<int, int>
#endif
#ifndef LD
#define LD long double
#endif


class KRTPushRelabelAdd {

    std::unordered_map<PII, int, boost::hash<PII>> g, cap;
    std::vector<int> d, excess, excess_hidden;
    std::unordered_set<PII, boost::hash<PII>> E_star;
    std::unordered_set<int> positiveEx;

    int N, M, S, T;
    std::vector<std::vector<KRTEdge>> G;

    EdgeDesignator *edgeDesignator;

    std::vector<dyn_item *> dynamicTreeNodes;
    std::unordered_map<dyn_item *, int> dynamicTreeIds;
    std::vector<std::unordered_set<int>> dynamicTreeChildren;

    Logger logger;

    // Util
    static PII rev(PII const &x) {
        return std::make_pair(x.second, x.first);
    }

    // Visible excess
    void maintain_positive_excess_data(int node) {
        if (excess[node] - excess_hidden[node] <= 0) positiveEx.erase(node);
        if (excess[node] - excess_hidden[node] > 0) positiveEx.insert(node);
    }

    int get_positive_excess_node() {
        positiveEx.erase(S);
        positiveEx.erase(T);
        if (!positiveEx.empty()) {
            return *positiveEx.begin();
        }
        return -1;
    }

    int get_excess_visible(int v) {
        return std::max(0, excess[v] - excess_hidden[v]);
    }

    // Flow
    int get_flow(PII const &e) {
        auto firstParent = dynamicTreeIds[dyn_find_father(dynamicTreeNodes[e.first])];
        auto secondParent = dynamicTreeIds[dyn_find_father(dynamicTreeNodes[e.second])];

        int res;
        if (firstParent == e.second) {
            // e in Ef
            res = cap[e] - dyn_find_value(dynamicTreeNodes[e.first]);
        } else if (secondParent == e.first) {
            // rev(e) in Ef
            res = -(cap[rev(e)] - dyn_find_value(dynamicTreeNodes[e.second]));
        } else {
            // e in E*
            res = g[e];
        }
        logger.log("f", e, res);
        return res;
    }

    void set_flow(PII const &e, int const &c) {
        logger.log("setFlow", e, c);
        g[e] = c;
        g[rev(e)] = -c;
    }

    void saturate(PII e) {
        logger.log("saturate", e, cap[e]);
        set_flow(e, cap[e]);

        excess[e.first] -= cap[e];
        excess[e.second] += cap[e];
        maintain_positive_excess_data(e.first);
        maintain_positive_excess_data(e.second);

        // sat push => adversary edge kill
        edgeDesignator->response_adversary(e.first, d[e.first], e.second, d[e.first] - 1, false);
    }

    // General
    void generic_init() {
        logger.log("genericInit");
        // Init helper vars
        for (int v = 1; v <= N; ++v) {
            for (auto const &nei: G[v]) {
                cap[std::make_pair(v, nei.dest)] = nei.cost;
                excess_hidden[v] += nei.cost;
                maintain_positive_excess_data(v);
            }
            d[v] = 0;
        }

        // Initialize Ef (link-cut tree)
        dynamicTreeNodes.resize(N + 1);
        dynamicTreeChildren.resize(N + 1);
        for (int v = 0; v <= N; ++v) {
            dynamicTreeNodes[v] = new dyn_item();
            dynamicTreeIds[dynamicTreeNodes[v]] = v;
            dyn_make_tree(dynamicTreeNodes[v], 0);
        }

        // Push flow from S at start
        // Sequence of relabels instead of simple d[S] = N; is needed to enforce adversary responses in the game
        for (int i = 0; i < N; ++i) relabel(S);
        for (auto const &nei: G[S]) {
            auto e = std::make_pair(S, nei.dest);
            add_edge(e);
        }
    }

    std::vector<PII > get_edges_list() {
        logger.log("initEdgesList");
        std::unordered_map<PII, int, boost::hash<PII>> edgesCosts;
        std::unordered_set<PII, boost::hash<PII>> edges;
        for (int i = 1; i <= N; ++i) {
            for (auto const &e: G[i]) {
                if (i == S || e.dest == S) continue;
                auto e_pair = std::make_pair(std::min(i, e.dest), std::max(i, e.dest));
                edgesCosts[e_pair] += e.cost;
                edges.insert(e_pair);
            }
        }
        std::vector<PII > edgesR(edges.begin(), edges.end());
        sort(edgesR.begin(), edgesR.end(), [&](PII const &a, PII const &b) {
            return edgesCosts[a] < edgesCosts[b];
        });
        return edgesR;
    }

    // TreePush/Relabel/AddEdge
    int get_bottleneck(int v) {
        int p = 0, q = std::numeric_limits<int>::max() - 1;
        dyn_item *root = dyn_find_root(dynamicTreeNodes[v]);
        while (p < q) {
            int mid = (p + q) / 2;
            dyn_item *bot = dyn_find_bottleneck(dynamicTreeNodes[v], mid);
            if (bot == root) {
                p = mid + 1;
            } else {
                q = mid;
            }
        }
        return p;
    }

    void link(PII e) {
        auto resCap = cap[e] - get_flow(e);
        logger.log("link", e, resCap);

        // Add edge to the tree
        dyn_link(dynamicTreeNodes[e.first], dynamicTreeNodes[e.second], resCap);
        dynamicTreeChildren[e.second].insert(e.first);
    }

    void cut(PII e) {
        logger.log("cut", e, -1);
        // Extract f value based on tree data
        auto resCap = dyn_find_value(dynamicTreeNodes[e.first]);
        set_flow(e, cap[e] - resCap);

        // Remove edge from the tree
        dyn_cut(dynamicTreeNodes[e.first]);
        dynamicTreeChildren[e.second].erase(e.first);

        // Sat push => adversary edge kill
        edgeDesignator->response_adversary(e.first, d[e.first], e.second, d[e.second], false);
    }

    void tree_push(int v, int vC) {
        logger.log("treePush", v);
        // Add current edge to the tree if needed
        if (dynamicTreeIds[dyn_find_root(dynamicTreeNodes[v])] == v) {
            link(std::make_pair(v, vC));
        }

        // Push flow
        auto pathMinimumResCap = get_bottleneck(v);
        auto excessVisible = get_excess_visible(v);
        auto delta = std::min(pathMinimumResCap, excessVisible);

        dyn_add_value(dynamicTreeNodes[v], -delta);

        // Update excess
        auto pathEndNode = dynamicTreeIds[dyn_find_root(dynamicTreeNodes[v])];
        excess[v] -= delta;
        excess[pathEndNode] += delta;
        maintain_positive_excess_data(v);
        maintain_positive_excess_data(pathEndNode);

        // Remove saturated edges
        dyn_item *start = dynamicTreeNodes[v];
        dyn_item *start_root = dyn_find_root(start);
        dyn_item *bot = dyn_find_bottleneck(start, 0);
        while (start != start_root && bot != start_root) {
            dyn_item *next_start = dyn_find_father(bot);
            cut(std::make_pair(dynamicTreeIds[bot], dynamicTreeIds[next_start]));
            start = next_start;
            start_root = dyn_find_root(start);
            bot = dyn_find_bottleneck(start, 0);
        }
    }

    void relabel(int v) {
        logger.log("relabel", v, d[v] + 1);
        auto childCopy = dynamicTreeChildren[v];
        for (auto const &c: childCopy) {
            cut(std::make_pair(c, v));
        }
        d[v]++;

        // relabel => remove node <v,d[v]> in the game
        edgeDesignator->response_adversary(-1, -1, v, d[v] - 1, true);

        // remove ineligible edges incident to v in the game
        for (auto const &w: G[v]) {
            auto e = std::make_pair(v, w.dest);
            // eligibility constraint
            if (!(cap[e] - get_flow(e) > 0 && d[v] == d[w.dest] + 1 && E_star.count(e))) {
                edgeDesignator->response_adversary(v, d[v], w.dest, d[v] - 1, false);
            }
        }
    }

    void add_edge(PII e) {
        logger.log("addEdge", e, -1);
        int v = e.first, w = e.second;

        // Add e to E*
        E_star.insert(e);
        E_star.insert(rev(e));

        // Subtract hidden excess
        excess_hidden[v] -= cap[e];
        excess_hidden[w] -= cap[rev(e)];
        maintain_positive_excess_data(v);
        maintain_positive_excess_data(w);

        // Saturate appropriate edge
        if (d[v] > d[w]) saturate(std::make_pair(v, w));
        else if (d[w] > d[v]) saturate(std::make_pair(w, v));
    }

public:
    // Initialize
    KRTPushRelabelAdd(int n, int m, int s, int t, const std::vector<std::vector<KRTEdge>> &g) : N(n), M(m), S(s), T(t),
                                                                                                G(g),
                                                                                                logger(Logger("KRT")) {
        d = excess = excess_hidden = std::vector<int>(n + 1);
        this->g.reserve(4 * m);
        cap.reserve(4 * m);
        E_star.rehash(4 * m);
        positiveEx.reserve(2 * n);
        edgeDesignator = new KRTEdgeDesignator(n, g);
    }

    // Run
    int run() {
        logger.log("run");

        generic_init();

        auto L = get_edges_list();
        while (!L.empty()) {
            auto e = L.back();
            L.pop_back();
            add_edge(e);

            int v;
            while ((v = get_positive_excess_node()) > 0) {
                auto currentEdge = edgeDesignator->ce(v, d[v]);

                logger.log("currentEdge", v, currentEdge);
                if (currentEdge == -1) {
                    relabel(v);
                } else {
                    tree_push(v, currentEdge);
                }
            }
        }
        return get_excess_visible(T);
    }
};

// API
int king_rao_tarjan(KRTGraph const &graph) {
    auto solver = KRTPushRelabelAdd(graph.getN(), graph.getM(), graph.getS(), graph.getT(), graph.getEdges());
    return solver.run();
}
