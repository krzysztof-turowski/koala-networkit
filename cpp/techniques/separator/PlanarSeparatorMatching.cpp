#include "techniques/separator/PlanarSeparatorMatching.hpp"
#include "matching/MaximumMatching.hpp"
#include "networkit/Globals.hpp"
#include "networkit/graph/EdgeUtils.hpp"
#include "techniques/separator/GraphUtils.hpp"
#include "techniques/separator/MISP.hpp"
#include <algorithm>
#include <queue>
#include <unordered_set>

namespace {
struct Action {
    int deg;
    NetworKit::node v = 0, u = 0, w = 0, x = 0;
};

std::vector<NetworKit::Edge> getExactMatching(NetworKit::Graph &graph) {
    Koala::EdmondsMaximumMatching m(graph, false);
    m.run();
    const auto &matchingMap = m.getMatching();
    std::vector<NetworKit::Edge> matchingVec;
    matchingVec.reserve(matchingMap.size());
    for (const auto &[u, v] : matchingMap) {
        if (u < v) {
            matchingVec.push_back({u, v});
        }
    }
    return matchingVec;
}

inline void insertEdgeIntoMatching(std::unordered_set<NetworKit::node> &matched_nodes,
                                   std::unordered_set<NetworKit::Edge> &S, NetworKit::Edge e) {
    S.insert(e);
    matched_nodes.insert(e.u);
    matched_nodes.insert(e.v);
}

inline void removeVertex(int &currentGraphSize, NetworKit::node v, std::vector<Action> &stk) {
    stk.push_back({0, v});
    currentGraphSize--;
}

inline void removeEdge(int &currentGraphSize, NetworKit::node v, std::queue<NetworKit::node> &Q,
                       std::vector<std::unordered_set<NetworKit::node>> &adj,
                       std::vector<Action> &stk) {
    NetworKit::node u = *adj[v].begin();
    stk.push_back({1, v, u});
    for (auto x : adj[u]) {
        adj[x].erase(u);
        if (adj[x].size() <= 2)
            Q.push(x);
    }
    adj[v].clear();
    adj[u].clear();
    currentGraphSize -= 2;
}

inline void contract(int &currentGraphSize, NetworKit::node v, std::queue<NetworKit::node> &Q,
                     std::vector<std::unordered_set<NetworKit::node>> &adj,
                     std::vector<Action> &stk) {
    auto it = adj[v].begin();
    NetworKit::node u = *it;
    NetworKit::node w = *std::next(it);
    NetworKit::node x = adj.size();
    adj.push_back({});
    std::unordered_set<NetworKit::node> seen;
    for (auto t : adj[u]) {
        if (t != v) {
            adj[t].erase(u);
            adj[t].insert(x);
            if (seen.find(t) == seen.end()) {
                adj[x].insert(t);
                seen.insert(t);
            }
        }
    }

    for (auto t : adj[w]) {
        if (t != v) {
            adj[t].erase(w);
            adj[t].insert(x);
            if (seen.find(t) == seen.end()) {
                adj[x].insert(t);
                seen.insert(t);
            }
        }
    }
    adj[v].clear();
    adj[u].clear();
    adj[w].clear();
    if (adj[x].size() <= 2)
        Q.push(x);

    stk.push_back({2, v, u, w, x});
    currentGraphSize -= 2;
}
} // namespace

namespace Koala {
PlanarSeparatorMatching::PlanarSeparatorMatching(NetworKit::Graph &G) : graph(G) {}

void PlanarSeparatorMatching::run() {
    reduce_procedure(graph);
    hasRun = true;
}

std::vector<NetworKit::Edge> PlanarSeparatorMatching::reduce_procedure(NetworKit::Graph &graph) {
    NetworKit::count n = graph.numberOfNodes();
    double loglog = std::max(1.0, std::log2(std::log2((double)std::max<size_t>(4, n))));

    int currentGraphSize = n;
    std::vector<Action> stk;
    std::queue<NetworKit::node> Q;
    std::vector<std::unordered_set<NetworKit::node>> adj(n);

    graph.forNodes([&](NetworKit::node v) {
        graph.forNeighborsOf(v, [&](NetworKit::node u) { adj[v].insert(u); });
    });

    graph.forNodes([&](NetworKit::node v) {
        if (adj[v].size() <= 2)
            Q.push(v);
    });

    std::unordered_set<NetworKit::Edge> S;
    std::unordered_set<NetworKit::node> matched_nodes;

    bool stop = false;
    while (!stop) {
        if (currentGraphSize <= loglog) {
            NetworKit::Graph induced = getInducedSubgraphFromAdj(adj);
            auto matching = getExactMatching(induced);
            for (auto e : matching) {
                insertEdgeIntoMatching(matched_nodes, S, e);
            }
            stop = true;
        } else if (!Q.empty()) {
            auto v = Q.front();
            Q.pop();

            if (adj[v].size() > 2)
                continue;

            if (adj[v].size() == 0) {
                removeVertex(currentGraphSize, v, stk);
            } else if (adj[v].size() == 1) {
                removeEdge(currentGraphSize, v, Q, adj, stk);
            } else if (adj[v].size() == 2) {
                contract(currentGraphSize, v, Q, adj, stk);
            }
        } else {
            // if there are no nodes with deg <= 2, we can use MISP algo
            NetworKit::Graph induced = getInducedSubgraphFromAdj(adj);

            double epsilon = loglog / std::max(currentGraphSize, 1);
            MISP<NetworKit::Edge> mispAlgo(
                induced, epsilon, [&](const NetworKit::Graph &c) -> std::vector<NetworKit::Edge> {
                    auto copyC = c;
                    auto matching = getExactMatching(copyC);
                    return matching;
                });
            mispAlgo.run();
            auto matching = mispAlgo.maximum_independent_set;
            for (auto e : matching) {
                insertEdgeIntoMatching(matched_nodes, S, e);
            }

            stop = true;
        }
    }

    while (!stk.empty()) {
        auto action = stk.back();
        stk.pop_back();

        if (action.deg == 0) {
            // just skip
        } else if (action.deg == 1) {
            if (matched_nodes.count(action.v) == 0 && matched_nodes.count(action.u) == 0) {
                insertEdgeIntoMatching(matched_nodes, S, NetworKit::Edge(action.v, action.u));
            }
        } else if (action.deg == 2) {
            // action.v -> degree 2 vertex matched to u and w
            // action.x -> virtual node combining u and w
            if (matched_nodes.count(action.x) == 0) {
                // Virtual node x was not matched, which means both u and w are free.
                // We can match v with u (or w). Let's match v with u.
                if (matched_nodes.count(action.v) == 0 && matched_nodes.count(action.u) == 0) {
                    insertEdgeIntoMatching(matched_nodes, S, NetworKit::Edge(action.v, action.u));
                }
            } else {
                // x was matched with some node t. Let's find t.
                NetworKit::node t = NetworKit::none;
                NetworKit::Edge matched_virtual_edge(NetworKit::none, NetworKit::none);
                for (const auto &e : S) {
                    if (e.u == action.x) {
                        t = e.v;
                        matched_virtual_edge = e;
                        break;
                    } else if (e.v == action.x) {
                        t = e.u;
                        matched_virtual_edge = e;
                        break;
                    }
                }

                if (t != NetworKit::none) {
                    // Remove the virtual edge from S and matched_nodes
                    S.erase(matched_virtual_edge);
                    matched_nodes.erase(action.x);

                    // Check if u is the one actually connected to t in G
                    if (graph.hasEdge(action.u, t)) {
                        // u matches with t; w is free to match with v
                        insertEdgeIntoMatching(matched_nodes, S, NetworKit::Edge(action.u, t));
                        if (matched_nodes.count(action.v) == 0 &&
                            matched_nodes.count(action.w) == 0) {
                            insertEdgeIntoMatching(matched_nodes, S,
                                                   NetworKit::Edge(action.v, action.w));
                        }
                    } else {
                        // w matches with t; u is free to match with v
                        insertEdgeIntoMatching(matched_nodes, S, NetworKit::Edge(action.w, t));
                        if (matched_nodes.count(action.v) == 0 &&
                            matched_nodes.count(action.u) == 0) {
                            insertEdgeIntoMatching(matched_nodes, S,
                                                   NetworKit::Edge(action.v, action.u));
                        }
                    }
                }
            }
        }
    }
    matching_set = std::vector<NetworKit::Edge>(S.begin(), S.end());
    hasRun = true;
}

} // namespace Koala
