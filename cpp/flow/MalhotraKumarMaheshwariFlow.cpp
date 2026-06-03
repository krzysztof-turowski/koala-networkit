#include <algorithm>
#include <limits>
#include <map>
#include <queue>
#include <vector>

#include "flow/MalhotraKumarMaheshwariFlow.hpp"
#include "graph/GraphTools.hpp"

namespace Koala {

const int UNREACHABLE = -2;

NetworKit::Edge MalhotraKumarMaheshwariFlow::reverse(const NetworKit::Edge &p) {
    return NetworKit::Edge(p.v, p.u);
}

void MalhotraKumarMaheshwariFlow::initialize() {
    graph.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        auto p = NetworKit::Edge(u, v);
        if (graph.addEdge(v, u, 0, true)) {
            capacity[reverse(p)] = 0;
        }
        flow[p] = 0;
        flow[reverse(p)] = 0;
        capacity[p] = w;
    });
}

bool MalhotraKumarMaheshwariFlow::build_level_graph() {
    graph.forNodes([&](NetworKit::node v) {
        level[v] = UNREACHABLE;
    });
    std::queue<NetworKit::node> q;
    q.push(source);
    level[source] = 0;
    while (!q.empty()) {
        NetworKit::node u = q.front();
        q.pop();

        graph.forNeighborsOf(u, [&](NetworKit::node w) {
            auto e = NetworKit::Edge(u, w);
            if (level[w] == UNREACHABLE && flow[e] < capacity[e]) {
                level[w] = level[u] + 1;
                q.push(w);
            }
        });
    }
    graph_stage = NetworKit::Graph(graph);
    return level[target] != UNREACHABLE;
}

void MalhotraKumarMaheshwariFlow::compute_potential() {
    graph.forNodes([&](NetworKit::node v) {
        in_potential[v] = 0;
        out_potential[v] = 0;
    });

    graph.forEdges([&](NetworKit::node x, NetworKit::node y) {
        auto e = NetworKit::Edge(x, y);
        if (level[e.u] + 1 == level[e.v]) {
            if (capacity[e] > flow[e]) {
                out_potential[e.u] += capacity[e] - flow[e];
                in_potential[e.v] += capacity[e] - flow[e];
            }
        }
    });

    in_potential[source] = std::numeric_limits<int>::max();
    out_potential[target] = std::numeric_limits<int>::max();
}

void MalhotraKumarMaheshwariFlow::push_forward(NetworKit::node u, NetworKit::edgeweight f) {
    if (u == target) {
        return;
    }
    std::queue<NetworKit::node> q;
    std::map<NetworKit::node, int> to_push;

    graph_stage.forNodes([&](NetworKit::node v) {
        to_push[v] = 0;
    });

    to_push[u] += f;
    q.push(u);

    while (!q.empty()) {
        NetworKit::node v = q.front();
        q.pop();
        if (to_push[v] == 0) {
            continue;
        }
        std::vector<NetworKit::Edge> edges_to_remove;
        graph_stage.forEdgesOf(v, [&](NetworKit::node w) {
            auto e = NetworKit::Edge(v, w);
            if (level[v] + 1 != level[e.v]) {
                return;
            }
            NetworKit::edgeweight can_be_pushed = std::min(capacity[e] - flow[e], to_push[v]);
            if (can_be_pushed == 0) {
                return;
            }
            if (e.v != target) {
                q.push(e.v);
            }
            flow[e] += can_be_pushed;
            flow[reverse(e)] -= can_be_pushed;
            if (capacity[e] - flow[e] == 0) {
                edges_to_remove.push_back(e);
            }
            in_potential[e.v] -= can_be_pushed;
            out_potential[v] -= can_be_pushed;
            to_push[v] -= can_be_pushed;
            to_push[e.v] += can_be_pushed;
        });
        for (const auto &e : edges_to_remove) {
            graph_stage.removeEdge(e.u, e.v);
            graph_stage.removeEdge(e.v, e.u);
        }
    }
}

void MalhotraKumarMaheshwariFlow::push_backward(NetworKit::node u, NetworKit::edgeweight f) {
    if (u == source) {
        return;
    }
    std::queue<NetworKit::node> q;
    std::map<NetworKit::node, int> to_push;

    graph_stage.forNodes([&](NetworKit::node v) {
        to_push[v] = 0;
    });

    to_push[u] += f;
    q.push(u);

    while (!q.empty()) {
        NetworKit::node v = q.front();
        q.pop();
        if (to_push[v] == 0) {
            continue;
        }
        std::vector<NetworKit::Edge> edges_to_remove;
        graph_stage.forInEdgesOf(v, [&](NetworKit::node w) {
            auto e = NetworKit::Edge(w, v);
            if (level[v] - 1 != level[e.u]) {
                return;
            }
            NetworKit::edgeweight can_be_pushed = std::min(capacity[e] - flow[e], to_push[v]);
            if (can_be_pushed == 0) {
                return;
            }
            if (e.u != source) {
                q.push(e.u);
            }
            flow[e] += can_be_pushed;
            flow[reverse(e)] -= can_be_pushed;
            if (capacity[e] - flow[e] == 0) {
                edges_to_remove.push_back(e);
            }
            out_potential[e.u] -= can_be_pushed;
            in_potential[v] -= can_be_pushed;
            to_push[v] -= can_be_pushed;
            to_push[e.u] += can_be_pushed;
        });
        for (const auto &e : edges_to_remove) {
            graph_stage.removeEdge(e.u, e.v);
            graph_stage.removeEdge(e.v, e.u);
        }
    }
}

void MalhotraKumarMaheshwariFlow::delete_node(NetworKit::node v) {
    graph_stage.forInEdgesOf(v, [&](NetworKit::node w) {
        auto e = NetworKit::Edge(w, v);
        if (level[v] - 1 == level[e.u]) {
            if (capacity[e] > flow[e]) {
                out_potential[e.u] -= capacity[e] - flow[e];
            }
        }
    });
    graph_stage.forEdgesOf(v, [&](NetworKit::node w) {
        auto e = NetworKit::Edge(v, w);
        if (level[v] + 1 == level[e.v]) {
            if (capacity[e] > flow[e]) {
                in_potential[e.v] -= capacity[e] - flow[e];
            }
        }
    });
    level[v] = UNREACHABLE;
    graph_stage.removeNode(v);
}

void MalhotraKumarMaheshwariFlow::run() {
    GraphTools::ensureDirectedGraph(graph);
    initialize();
    int total_flow = 0;
    while (build_level_graph()) {
        compute_potential();
        while (true) {
            NetworKit::node u;
            int minimum = std::numeric_limits<int>::max();

            graph_stage.forNodes([&](NetworKit::node v) {
                if (level[v] == UNREACHABLE) {
                    return;
                }
                int pot = std::min(out_potential[v], in_potential[v]);
                if (pot < minimum) {
                    u = v;
                    minimum = pot;
                }
            });
            if (minimum == 0) {
                delete_node(u);
                continue;
            }
            if (minimum == std::numeric_limits<int>::max()) {
                break;
            }

            push_forward(u, minimum);
            push_backward(u, minimum);
            total_flow += minimum;
            delete_node(u);
        }
    }
    flow_size = total_flow;
    hasRun = true;
}
}  // namespace Koala
