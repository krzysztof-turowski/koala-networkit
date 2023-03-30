#include <bits/stdc++.h>

#include <boost/functional/hash.hpp>

#include <networkit/graph/GraphTools.hpp>

#include "krt_edge_designator.hpp"
#include "logger.hpp"

#include "dynamic_trees/dyn_tree.h"

#ifndef PII
#define PII std::pair<int, int>
#endif
#ifndef LD
#define LD long double
#endif


using edge = std::pair<NetworKit::node, NetworKit::node>;

class KRTPushRelabelAdd {
    const std::optional<NetworKit::Graph> graph;

    std::unordered_map<PII, int, boost::hash<PII>> flow, capacity;
    std::vector<int> d, excess, excess_hidden;
    std::unordered_set<PII, boost::hash<PII>> E_star;
    std::unordered_set<int> positive_excess;

    int S, T;

    KRTEdgeDesignator *edge_designator;

    std::vector<dyn_item *> dynamicTreeNodes;
    std::unordered_map<dyn_item *, int> dynamicTreeIds;
    std::vector<std::unordered_set<int>> dynamicTreeChildren;

    Logger logger;

edge reverse(const edge &p) {
    return std::make_pair(p.second, p.first);
}

int get_visible_excess(NetworKit::node v) {
    return std::max(0, excess[v] - excess_hidden[v]);
}

void update_positive_excess_map(NetworKit::node v) {
    if (get_visible_excess(v) > 0) {
        positive_excess.insert(v);
    } else {
        positive_excess.erase(v);
    }
}

    int get_positive_excess_node() {
        positive_excess.erase(S);
        positive_excess.erase(T);
        if (!positive_excess.empty()) {
            return *positive_excess.begin();
        }
        return -1;
    }

    // Flow
    int get_flow(PII const &e) {
        auto firstParent = dynamicTreeIds[dyn_find_father(dynamicTreeNodes[e.first])];
        auto secondParent = dynamicTreeIds[dyn_find_father(dynamicTreeNodes[e.second])];

        int res;
        if (firstParent == e.second) {
            // e in Ef
            res = capacity[e] - dyn_find_value(dynamicTreeNodes[e.first]);
        } else if (secondParent == e.first) {
            // reverse(e) in Ef
            res = -(capacity[reverse(e)] - dyn_find_value(dynamicTreeNodes[e.second]));
        } else {
            // e in E*
            res = flow[e];
        }
        logger.log("f", e, res);
        return res;
    }

void set_flow(const edge &e, int c) {
    logger.log("setFlow", e, c);
    flow[e] = c, flow[reverse(e)] = -c;
}

void saturate(const edge &e) {
    logger.log("saturate", e, capacity[e]);
    set_flow(e, capacity[e]);
    excess[e.first] -= capacity[e], excess[e.second] += capacity[e];
    update_positive_excess_map(e.first);
    update_positive_excess_map(e.second);

    // sat push => adversary edge kill
    edge_designator->response_adversary(e.first, d[e.first], e.second, d[e.first] - 1, false);
}

void initialize() {
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        capacity[std::make_pair(u + 1, v + 1)] = w;
        excess_hidden[u + 1] += w;
        update_positive_excess_map(u + 1);
    });
    graph->forNodes([&](NetworKit::node v) {
        d[v + 1] = 0;
    });

    // dynamic_tree.initialize();
    // Initialize Ef (link-cut tree)
    dynamicTreeNodes.resize(graph->numberOfNodes() + 1);
    dynamicTreeChildren.resize(graph->numberOfNodes() + 1);
    for (int v = 0; v <= graph->numberOfNodes(); ++v) {
        dynamicTreeNodes[v] = new dyn_item();
        dynamicTreeIds[dynamicTreeNodes[v]] = v;
        dyn_make_tree(dynamicTreeNodes[v], 0);
    }

    for (int i = 0; i < graph->numberOfNodes(); i++) {
        relabel(S);
    }
    graph->forNeighborsOf(S - 1, [&](NetworKit::node v) {
        add_edge(std::make_pair(S, v + 1));
    });
}

    std::vector<PII > get_edges_list() {
        logger.log("initEdgesList");
        std::unordered_map<PII, int, boost::hash<PII>> edge_costs;
        std::unordered_set<PII, boost::hash<PII>> edges;
        graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
            if (u + 1 == S || v + 1 == S) {
                return;
            }
            auto e = std::make_pair(std::min(u + 1, v + 1), std::max(u + 1, v + 1));
            edge_costs[e] += w;
            edges.insert(e);
        });
        std::vector<PII> sorted_edges(edges.begin(), edges.end());
        sort(sorted_edges.begin(), sorted_edges.end(), [&](PII const &a, PII const &b) {
            return edge_costs[a] < edge_costs[b];
        });
        return sorted_edges;
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
        auto resCap = capacity[e] - get_flow(e);
        logger.log("link", e, resCap);

        // Add edge to the tree
        dyn_link(dynamicTreeNodes[e.first], dynamicTreeNodes[e.second], resCap);
        dynamicTreeChildren[e.second].insert(e.first);
    }

    void cut(PII e) {
        logger.log("cut", e, -1);
        // Extract f value based on tree data
        auto resCap = dyn_find_value(dynamicTreeNodes[e.first]);
        set_flow(e, capacity[e] - resCap);

        // Remove edge from the tree
        dyn_cut(dynamicTreeNodes[e.first]);
        dynamicTreeChildren[e.second].erase(e.first);

        // Sat push => adversary edge kill
        edge_designator->response_adversary(e.first, d[e.first], e.second, d[e.second], false);
    }

    void tree_push(int v, int vC) {
        logger.log("treePush", v);
        // Add current edge to the tree if needed
        if (dynamicTreeIds[dyn_find_root(dynamicTreeNodes[v])] == v) {
            link(std::make_pair(v, vC));
        }

        // Push flow
        auto pathMinimumResCap = get_bottleneck(v);
        auto excessVisible = get_visible_excess(v);
        auto delta = std::min(pathMinimumResCap, excessVisible);

        dyn_add_value(dynamicTreeNodes[v], -delta);

        // Update excess
        auto pathEndNode = dynamicTreeIds[dyn_find_root(dynamicTreeNodes[v])];
        excess[v] -= delta;
        excess[pathEndNode] += delta;
        update_positive_excess_map(v);
        update_positive_excess_map(pathEndNode);

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
        edge_designator->response_adversary(-1, -1, v, d[v] - 1, true);

        // remove ineligible edges incident to v in the game
        graph->forNeighborsOf(v - 1, [&](NetworKit::node w) {
            auto e = std::make_pair(v, w + 1);
            // eligibility constraint
            if (!(capacity[e] - get_flow(e) > 0 && d[v] == d[w + 1] + 1 && E_star.count(e))) {
                edge_designator->response_adversary(v, d[v], w + 1, d[v] - 1, false);
            }
        });
    }

void add_edge(const edge &e) {
    const edge &e_rev = reverse(e);
    E_star.insert(e), E_star.insert(e_rev);
    excess_hidden[e.first] -= capacity[e], excess_hidden[e.second] -= capacity[e_rev];
    update_positive_excess_map(e.first);
    update_positive_excess_map(e.second);

    int &d_first = d[e.first], &d_second = d[e.second];
    if (d_first > d_second) {
        saturate(e);
    } else if (d_second > d_first) {
        saturate(e_rev);
    }
}

public:
    // Initialize
    KRTPushRelabelAdd(
            const std::optional<NetworKit::Graph> &graph,
            NetworKit::node source, NetworKit::node target) : graph(graph), S(source + 1), T(target + 1), logger(Logger("KRT")) {
        d = excess = excess_hidden = std::vector<int>(graph->numberOfNodes() + 1);
        flow.reserve(4 * graph->numberOfEdges());
        capacity.reserve(4 * graph->numberOfEdges());
        E_star.rehash(4 * graph->numberOfEdges());
        positive_excess.reserve(2 * graph->numberOfNodes());
        edge_designator = new KRTEdgeDesignator(graph->numberOfNodes(), graph);
    }

    // Run
    int run() {
        logger.log("run");
        initialize();
        auto L = get_edges_list();
        while (!L.empty()) {
            auto e = L.back();
            L.pop_back();
            add_edge(e);

            int v;
            while ((v = get_positive_excess_node()) > 0) {
                auto currentEdge = edge_designator->ce(v, d[v]);

                logger.log("currentEdge", v, currentEdge);
                if (currentEdge == -1) {
                    relabel(v);
                } else {
                    tree_push(v, currentEdge);
                }
            }
        }
        return get_visible_excess(T);
    }
};

// API
int king_rao_tarjan(
      const std::optional<NetworKit::Graph> &graph,
      NetworKit::node source, NetworKit::node target) {
    auto solver = KRTPushRelabelAdd(graph, source, target);
    return solver.run();
}
