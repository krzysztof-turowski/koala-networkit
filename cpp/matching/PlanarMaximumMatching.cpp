#include <algorithm>
#include <queue>
#include <unordered_set>
#include <vector>
#include "matching/MaximumMatching.hpp"
#include "networkit/Globals.hpp"
#include "networkit/graph/EdgeUtils.hpp"
#include "techniques/separator/MISP.hpp"
#include <boost/parameter/aux_/pp_impl/match.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace {} // namespace

namespace Koala {
PlanarSeparatorMatching::PlanarSeparatorMatching(NetworKit::Graph &G)
    : MaximumCardinalityMatching(G) {}

void PlanarSeparatorMatching::run() {
    reduce_procedure(graph);
    hasRun = true;
}

std::vector<NetworKit::Edge>
PlanarSeparatorMatching::get_exact_matching(const NetworKit::Graph &G) {
    NetworKit::Graph copyG = G;
    Koala::EdmondsMaximumMatching m(copyG, false);
    m.run();
    std::vector<NetworKit::Edge> matchingVec;
    for (auto [u, v] : m.getMatching()) {
        if (u == NetworKit::none || v == NetworKit::none)
            continue;

        if (u < v) {
            matchingVec.push_back(NetworKit::Edge(u, v));
        }
    }
    return matchingVec;
}

void PlanarSeparatorMatching::contract(int &currentGraphSize, NetworKit::node v,
                                       std::queue<NetworKit::node> &Q, NetworKit::Graph &G,
                                       std::vector<Action> &stk) {
    NetworKit::node u = NetworKit::none;
    NetworKit::node w = NetworKit::none;
    // v has degree 2, so the first neighbor is u and second is w
    G.forNeighborsOf(v, [&](NetworKit::node t) {
        if (u == NetworKit::none) {
            u = t;
        } else {
            w = t;
        }
    });

    NetworKit::node larger = (G.degree(u) >= G.degree(w)) ? u : w;
    NetworKit::node smaller = (larger == u) ? w : u;
    Action action;
    action.type = ReductionType::DEG_2;
    action.v = v;
    action.larger = larger;
    action.smaller = smaller;

    std::unordered_set<NetworKit::node> neighborsOfLarger;
    G.forNeighborsOf(larger, [&](NetworKit::node t) { neighborsOfLarger.insert(t); });

    G.forNeighborsOf(smaller, [&](NetworKit::node t) {
        if (t != v && t != larger) {
            if (!neighborsOfLarger.contains(t)) {
                G.addEdge(larger, t);
                neighborsOfLarger.insert(t);
                action.moved_nodes.push_back(t);
            }
        }
    });

    if (G.degree(larger) <= 2) {
        Q.push(larger);
    }

    G.removeNode(v);
    G.removeNode(smaller);
    stk.push_back(action);
    currentGraphSize -= 2;
}

void PlanarSeparatorMatching::remove_edge(int &currentGraphSize, NetworKit::node v,
                                          std::queue<NetworKit::node> &Q, NetworKit::Graph &G,
                                          std::vector<Action> &stk) {
    NetworKit::node u = NetworKit::none;
    G.forNeighborsOf(v, [&](NetworKit::node t) { u = t; });

    std::vector<NetworKit::node> neighbors_of_u;
    G.forNeighborsOf(u, [&](NetworKit::node t) {
        if (t != v) {
            neighbors_of_u.push_back(t);
        }
    });

    G.removeNode(u);
    G.removeNode(v);

    for (auto t : neighbors_of_u) {
        if (G.degree(t) <= 2)
            Q.push(t);
    }

    Action action;
    action.type = ReductionType::DEG_1;
    action.v = v;
    action.u = u;
    stk.push_back(action);
    currentGraphSize -= 2;
}

void PlanarSeparatorMatching::remove_vertex(int &currentGraphSize, NetworKit::node v,
                                            NetworKit::Graph &G, std::vector<Action> &stk) {
    G.removeNode(v);
    Action action;
    action.type = ReductionType::DEG_0;
    action.v = v;
    stk.push_back(action);
    currentGraphSize--;
}

void PlanarSeparatorMatching::insert_edge_into_matching(
    std::unordered_set<NetworKit::node> &matched_nodes, std::unordered_set<NetworKit::Edge> &S,
    NetworKit::Edge e) {
    S.insert(e);
    matched_nodes.insert(e.u);
    matched_nodes.insert(e.v);
}

void PlanarSeparatorMatching::reduce_procedure(NetworKit::Graph &graph) {
    NetworKit::Graph localGraph = graph;
    NetworKit::count n = localGraph.numberOfNodes();

    double loglog = std::max(1.0, std::log2(std::log2((double)std::max<size_t>(4, n))));

    int currentGraphSize = n;
    std::vector<Action> stk;
    std::queue<NetworKit::node> Q;
    localGraph.forNodes([&](NetworKit::node v) {
        if (localGraph.degree(v) <= 2)
            Q.push(v);
    });

    std::unordered_set<NetworKit::Edge> S;
    std::unordered_set<NetworKit::node> matched_nodes;

    bool stop = false;
    while (!stop) {
        if (currentGraphSize <= loglog) {
            auto localMatching = get_exact_matching(localGraph);
            if (!localGraph.hasEdgeIds()) {
                localGraph.indexEdges();
            }
            for (auto e : localMatching) {
                insert_edge_into_matching(matched_nodes, S, NetworKit::Edge(e.u, e.v));
            }
            stop = true;
        } else if (!Q.empty()) {
            auto v = Q.front();
            Q.pop();

            if (!localGraph.hasNode(v)) {
                continue;
            }
            int deg_v = localGraph.degree(v);

            if (deg_v > 2)
                continue;

            if (deg_v == 0) {
                remove_vertex(currentGraphSize, v, localGraph, stk);
            } else if (deg_v == 1) {
                remove_edge(currentGraphSize, v, Q, localGraph, stk);
            } else if (deg_v == 2) {
                contract(currentGraphSize, v, Q, localGraph, stk);
            }
        } else {
            double epsilon = loglog / std::max(currentGraphSize, 1);
            MISP<NetworKit::Edge> mispAlgo(
                localGraph, epsilon,
                [&](const NetworKit::Graph &c) -> std::vector<NetworKit::Edge> {
                    if (c.numberOfEdges() == 0)
                        return {};
                    NetworKit::Graph copyG = c;
                    if (!copyG.hasEdgeIds()) {
                        copyG.indexEdges();
                    }
                    auto localMatching = get_exact_matching(copyG);
                    return localMatching;
                });
            mispAlgo.run();
            auto localMatching = mispAlgo.getIndependentSet();
            for (auto e : localMatching) {
                insert_edge_into_matching(matched_nodes, S, NetworKit::Edge(e.u, e.v));
            }
            stop = true;
        }
    }
    matching.clear();
    graph.forNodes([&](NetworKit::node v) { matching[v] = NetworKit::none; });
    for (auto e : S) {
        matching[e.u] = e.v;
        matching[e.v] = e.u;
    }

    while (!stk.empty()) {
        auto action = stk.back();
        stk.pop_back();

        if (action.type == ReductionType::DEG_0) {
            // just skip
        } else if (action.type == ReductionType::DEG_1) {
            matching[action.u] = action.v;
            matching[action.v] = action.u;
        } else if (action.type == ReductionType::DEG_2) {
            auto v = action.v;
            auto larger = action.larger;
            auto smaller = action.smaller;

            NetworKit::node matchedWithLarger = NetworKit::none;
            if (matching.count(larger)) {
                matchedWithLarger = matching[larger];
            }

            bool edgeBelongedToSmaller = false;
            for (auto neighbor : action.moved_nodes) {
                if (neighbor == matchedWithLarger) {
                    edgeBelongedToSmaller = true;
                    break;
                }
            }

            if (edgeBelongedToSmaller) {
                matching.erase(larger);
                matching[smaller] = matchedWithLarger;
                matching[matchedWithLarger] = smaller;

                matching[larger] = v;
                matching[v] = larger;
            } else {
                matching[smaller] = v;
                matching[v] = smaller;
            }
        }
    }
}

} // namespace Koala
