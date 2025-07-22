#include <limits>
#include <algorithm>
#include <queue>
#include <flow/BoykovKolmogorovFlow.hpp>

using edge = NetworKit::Edge;
constexpr int SOURCE = 0;
constexpr int TARGET = 1;
constexpr int FREE   = 2;
namespace Koala {

edge BoykovKolmogorovFlow::reverse(const edge &p) {
    return NetworKit::Edge(p.v, p.u);
}

void BoykovKolmogorovFlow::initialize() {
    V = graph->numberOfNodes();
    graph->forNodes([&](NetworKit::node v) {
        tree[v] = FREE;
        parent[v] = NetworKit::none;
    });
    tree[source] = SOURCE;
    tree[target] = TARGET;
    active.push(source);
    active.push(target);
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

int  BoykovKolmogorovFlow::tree_capacity(NetworKit::node p, NetworKit::node q) {
    edge forward = NetworKit::Edge(p, q);
    edge backward = NetworKit::Edge(q, p);
    if (tree[p] == SOURCE) {
        return capacity[forward] - flow[forward];
    }
    if (tree[p] == TARGET) {
        return capacity[backward] - flow[backward];
    }
    return 0;
}

bool  BoykovKolmogorovFlow::grow() {
    while (!active.empty()) {
        NetworKit::node v = active.front();
        active.pop();
        bool foundPath = false;
        if (tree[v] == FREE) {
            continue;
        }
        graph->forNeighborsOf(v, [&](NetworKit::node w) {
            if (foundPath) return;
            edge e = NetworKit::Edge(v, w);
            if (tree_capacity(v, w) > 0) {
                if (tree[w] == FREE) {
                    tree[w] = tree[v];
                    parent[w] = v;
                    active.push(w);
                } else if (tree[w] != tree[v]) {
                    spath = tree[v] == SOURCE ? v : w;
                    tpath = tree[w] == TARGET ? w : v;
                    foundPath = true;
                    active.push(v);
                    return;
                }
            }
        });
        if (foundPath) return true;
    }
    return false;
}

int BoykovKolmogorovFlow::augment() {
    std::vector<edge> path;
    int bottleneck = std::numeric_limits<int>::max();

    NetworKit::node u = spath;
    while (u != source) {
        NetworKit::node p = parent[u];
        edge e = NetworKit::Edge(p, u);
        bottleneck = std::min(bottleneck, capacity[e] - flow[e]);
        path.push_back(e);
        u = p;
    }

    edge middle = NetworKit::Edge(spath, tpath);
    bottleneck = std::min(bottleneck, capacity[middle] - flow[middle]);
    path.push_back(middle);

    u = tpath;
    while (u != target) {
        NetworKit::node p = parent[u];
        edge e = NetworKit::Edge(u, p);  // reverse direction
        bottleneck = std::min(bottleneck, capacity[e] - flow[e]);
        path.push_back(e);
        u = p;
    }

    for (edge e : path) {
        flow[e] += bottleneck;
        flow[reverse(e)] -= bottleneck;
    }

    u = spath;
    while (u != source) {
        NetworKit::node p = parent[u];
        edge e = NetworKit::Edge(p, u);
        if (capacity[e] - flow[e] == 0) {
            if (tree[p] == SOURCE && tree[u] == SOURCE) {
                parent[u] = NetworKit::none;
                orphan.push(u);
            }
        }
        u = p;
    }

    u = tpath;
    while (u != target) {
        NetworKit::node p = parent[u];
        edge e = NetworKit::Edge(u, p);  // reverse direction
        if (capacity[e] - flow[e] == 0) {
            if (tree[p] == TARGET && tree[u] == TARGET) {
                parent[u] = NetworKit::none;
                orphan.push(u);
            }
        }
        u = p;
    }
    return bottleneck;
}

bool BoykovKolmogorovFlow::origin(NetworKit::node v) {
    NetworKit::node u = v;
    while (true) {
        if (parent[u] == source || parent[u] == target) {
            return true;
        }
        if (parent[u] == NetworKit::none) {
            return false;
        }
        u = parent[u];
    }
}

void BoykovKolmogorovFlow::adopt() {
    while (!orphan.empty()) {
        NetworKit::node p = orphan.front();
        orphan.pop();
        bool found_new_parent = false;

        graph->forNeighborsOf(p, [&](NetworKit::node q) {
            if (found_new_parent)return;

            if (tree[q] != tree[p]) return;
            if (tree_capacity(q, p) <= 0) return;
            if (!origin(q)) return;

            parent[p] = q;
            found_new_parent = true;
        });

        if (!found_new_parent) {
            graph->forNeighborsOf(p, [&](NetworKit::node q) {
                if (tree[q] != tree[p]) return;

                if (tree_capacity(q, p) > 0) {
                    active.push(q);
                }

                if (parent[q] == p) {
                    parent[q] = NetworKit::none;
                    orphan.push(q);
                }
            });
            tree[p] = FREE;
        }
    }
}

void BoykovKolmogorovFlow::run() {
    initialize();
    int totalflow = 0;
    int iteration = 0;

    while (true) {
        if (!grow()) {
            break;
        }

        int f = augment();
        totalflow += f;
        adopt();
    }

    flow_size = totalflow;
    hasRun = true;
}
}  // namespace Koala
