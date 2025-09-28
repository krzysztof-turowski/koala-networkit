#include <limits>
#include <algorithm>
#include <queue>
#include <flow/MalhotraKumarMaheshwariFlow.hpp>

using edge = NetworKit::Edge;

namespace Koala {

const int UNREACHABLE = -2;

edge MalhotraKumarMaheshwariFlow::reverse(const edge &p) {
    return NetworKit::Edge(p.v, p.u);
}

void MalhotraKumarMaheshwariFlow::initialize() {
    V = graph->numberOfNodes();
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        auto p = NetworKit::Edge(u, v);
        if (graph->addEdge(v, u, 0, true)) {
            capacity[reverse(p)] = 0;
        }
        flow[p] = 0;
        flow[reverse(p)] = 0;
        capacity[p] = w;
    });
}

bool MalhotraKumarMaheshwariFlow::buildLevelGraph() {
    graph->forNodes([&](NetworKit::node v) {
        level[v] = UNREACHABLE;
    });
    std::queue<NetworKit::node> q;
    q.push(source);
    level[source] = 0;
    while (!q.empty()) {
        NetworKit::node u = q.front();
        q.pop();

        graph->forNeighborsOf(u, [&](NetworKit::node w) {
            auto e = NetworKit::Edge(u, w);
            if (level[w] == UNREACHABLE && flow[e] < capacity[e]) {
                level[w] = level[u] + 1;
                q.push(w);
            }
        });
    }
    graph_stage = NetworKit::Graph(*graph);
    return level[target] != UNREACHABLE;
}

void MalhotraKumarMaheshwariFlow::computePotential() {
    graph->forNodes([&](NetworKit::node v) {
        inPotential[v] = 0;
        outPotential[v] = 0;
    });

    graph->forEdges([&](NetworKit::node x, NetworKit::node y, NetworKit::edgeweight weight) {
        edge e = NetworKit::Edge(x, y);
        if (level[e.u] + 1 == level[e.v]) {
            if (capacity[e] > flow[e]) {
                outPotential[e.u]  +=  capacity[e] - flow[e];
                inPotential[e.v]  +=  capacity[e] - flow[e];
            }
        }
    });

    inPotential[source] = std::numeric_limits<int>::max();
    outPotential[target] = std::numeric_limits<int>::max();
}

void MalhotraKumarMaheshwariFlow::pushForward(NetworKit::node u, NetworKit::edgeweight f) {
    if (u == target) {
        return;
    }
    std::queue<NetworKit::node> q;
    std::map<NetworKit::node, int> to_push;

    graph_stage.forNodes([&](NetworKit::node v) {
        to_push[v] = 0;
    });

    to_push[u]  +=  f;
    q.push(u);

    while (!q.empty()) {
        NetworKit::node v = q.front();
        q.pop();
        if (to_push[v] == 0) {
            continue;
        }
        std::vector<edge> edgesToRemove;
        graph_stage.forEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight) {
            edge e = NetworKit::Edge(v, w);
            if (level[v] + 1 != level[e.v]) {
                return;
            }
            NetworKit::edgeweight can_be_pushed = std::min(capacity[e]-flow[e], to_push[v]);
            if (can_be_pushed == 0) {
                return;
            }
            if (e.v != target) {
                q.push(e.v);
            }
            flow[e] += can_be_pushed;
            flow[reverse(e)] -= can_be_pushed;
            if (capacity[e] - flow[e] == 0) {
                edgesToRemove.push_back(e);
            }
            inPotential[e.v] -= can_be_pushed;
            outPotential[v] -= can_be_pushed;
            to_push[v] -= can_be_pushed;
            to_push[e.v] += can_be_pushed;
        });
        for (const edge& e : edgesToRemove) {
            graph_stage.removeEdge(e.u, e.v);
            graph_stage.removeEdge(e.v, e.u);
        }
    }
}

void MalhotraKumarMaheshwariFlow::pushBackward(NetworKit::node u, NetworKit::edgeweight f) {
    if (u == source) {
        return;
    }
    std::queue<NetworKit::node> q;
    std::map<NetworKit::node, int> to_push;

    graph_stage.forNodes([&](NetworKit::node v) {
        to_push[v] = 0;
    });

    to_push[u]  +=  f;
    q.push(u);

    while (!q.empty()) {
        NetworKit::node v = q.front();
        q.pop();
        if (to_push[v] == 0) {
            continue;
        }
        std::vector<edge> edgesToRemove;
        graph_stage.forInEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight) {
            edge e = NetworKit::Edge(w, v);
            if (level[v] - 1 != level[e.u]) {
                return;
            }
            NetworKit::edgeweight can_be_pushed = std::min(capacity[e]-flow[e], to_push[v]);
            if (can_be_pushed == 0) {
                return;
            }
            if (e.u != source) {
                q.push(e.u);
            }
            flow[e] += can_be_pushed;
            flow[reverse(e)] -= can_be_pushed;
            if (capacity[e] - flow[e] == 0) {
                edgesToRemove.push_back(e);
            }
            outPotential[e.u] -= can_be_pushed;
            inPotential[v] -= can_be_pushed;
            to_push[v] -= can_be_pushed;
            to_push[e.u] += can_be_pushed;
        });
        for (const edge& e : edgesToRemove) {
            graph_stage.removeEdge(e.u, e.v);
            graph_stage.removeEdge(e.v, e.u);
        }
    }
}

void MalhotraKumarMaheshwariFlow::deleteNode(NetworKit::node v) {
     graph_stage.forInEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight) {
        edge e = NetworKit::Edge(w, v);
        if (level[v] - 1 == level[e.u]) {
            if (capacity[e] > flow[e]) {
                outPotential[e.u]  -=  capacity[e] - flow[e];
            }
        }
    });
    graph_stage.forEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight) {
        edge e = NetworKit::Edge(v, w);
        if (level[v] + 1 == level[e.v]) {
            if (capacity[e] > flow[e]) {
                inPotential[e.v]  -=  capacity[e] -flow[e];
            }
        }
    });
    level[v] = UNREACHABLE;
    graph_stage.removeNode(v);
}

void MalhotraKumarMaheshwariFlow::run() {
    initialize();
    int totalFlow = 0;
    while (buildLevelGraph()) {
        computePotential();
        while (true) {
            NetworKit::node u;
            int minimum = std::numeric_limits<int>::max();

            graph_stage.forNodes([&](NetworKit::node v) {
                if (level[v] == UNREACHABLE) {
                    return;
                }
                int pot = std::min(outPotential[v], inPotential[v]);
                if (pot < minimum) {
                    u = v;
                    minimum = pot;
                }
            });
            if (minimum == 0) {
                deleteNode(u);
                continue;
            }
            if (minimum == std::numeric_limits<int>::max()) {
                break;
            }

            pushForward(u, minimum);
            pushBackward(u, minimum);
            totalFlow += minimum;
            deleteNode(u);
        }
    }
    flow_size = totalFlow;
    hasRun = true;
}
}  // namespace Koala
