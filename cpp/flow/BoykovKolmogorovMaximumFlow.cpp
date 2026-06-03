#include <algorithm>
#include <limits>
#include <queue>
#include <vector>

#include "flow/BoykovKolmogorovMaximumFlow.hpp"
#include "graph/GraphTools.hpp"

namespace Koala {

static NetworKit::Edge reverse(const NetworKit::Edge &p) {
    return NetworKit::Edge(p.v, p.u);
}

void BoykovKolmogorovMaximumFlow::initialize(std::queue<NetworKit::node> &active) {
    graph.forNodes([&](NetworKit::node v) {
        tree[v] = NodeType::FREE;
        parent[v] = NetworKit::none;
    });
    tree[source] = NodeType::SOURCE;
    tree[target] = NodeType::TARGET;
    active.push(source);
    active.push(target);
    graph.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        auto p = NetworKit::Edge(u, v);
        if (graph.addEdge(v, u, 0, true)) {
            capacity[reverse(p)] = 0;
        }
        flow[p] = 0;
        flow[reverse(p)] = 0;
        capacity[p] = w;
    });
    flow_size = 0;
}

int BoykovKolmogorovMaximumFlow::tree_capacity(NetworKit::node p, NetworKit::node q) {
    auto forward = NetworKit::Edge(p, q);
    auto backward = NetworKit::Edge(q, p);
    if (tree[p] == NodeType::SOURCE) {
        return capacity[forward] - flow[forward];
    }
    if (tree[p] == NodeType::TARGET) {
        return capacity[backward] - flow[backward];
    }
    return 0;
}

std::optional<NetworKit::Edge> BoykovKolmogorovMaximumFlow::grow(
    std::queue<NetworKit::node> &active) {
    while (!active.empty()) {
        NetworKit::node v = active.front();
        active.pop();
        if (tree[v] == NodeType::FREE) {
            continue;
        }

        const auto neighbors = graph.neighborRange(v);
        const auto found_vertex = std::find_if(
            neighbors.begin(), neighbors.end(),
            [&](NetworKit::node u) {
                if (tree_capacity(v, u) <= 0) {
                    return false;
                }
                if (tree[u] == NodeType::FREE) {
                    tree[u] = tree[v];
                    parent[u] = v;
                    active.push(u);
                    return false;
                }
                return tree[u] != tree[v];
            });

        if (found_vertex != neighbors.end()) {
            const NetworKit::node u = *found_vertex;
            const NetworKit::node spath = tree[v] == NodeType::SOURCE ? v : u;
            const NetworKit::node tpath = tree[u] == NodeType::TARGET ? u : v;
            active.push(v);
            return NetworKit::Edge(spath, tpath);
        }
    }
    return std::nullopt;
}

int BoykovKolmogorovMaximumFlow::augment(
    const NetworKit::Edge &middle, std::queue<NetworKit::node> &orphan) {
    std::vector<NetworKit::Edge> path;
    auto add_path = [&](NetworKit::node u, NetworKit::node root, bool reverse_edges) {
        while (u != root) {
            const NetworKit::node p = parent[u];
            path.emplace_back(reverse_edges ? NetworKit::Edge(u, p) : NetworKit::Edge(p, u));
            u = p;
        }
    };

    add_path(middle.u, source, false);
    path.push_back(middle);
    add_path(middle.v, target, true);

    int bottleneck = std::numeric_limits<int>::max();
    for (const auto &e : path) {
        bottleneck = std::min(bottleneck, capacity[e] - flow[e]);
    }

    for (const auto &e : path) {
        flow[e] += bottleneck;
        flow[reverse(e)] -= bottleneck;
        if (capacity[e] - flow[e] != 0) {
            continue;
        }
        if (tree[e.v] == NodeType::SOURCE && parent[e.v] == e.u) {
            parent[e.v] = NetworKit::none;
            orphan.push(e.v);
        }
        if (tree[e.u] == NodeType::TARGET && parent[e.u] == e.v) {
            parent[e.u] = NetworKit::none;
            orphan.push(e.u);
        }
    }
    return bottleneck;
}

bool BoykovKolmogorovMaximumFlow::origin(NetworKit::node v) {
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

void BoykovKolmogorovMaximumFlow::adopt(
    std::queue<NetworKit::node> &active, std::queue<NetworKit::node> &orphan) {
    while (!orphan.empty()) {
        auto p = orphan.front();
        orphan.pop();

        const auto neighbors = graph.neighborRange(p);
        const auto new_parent = std::find_if(
            neighbors.begin(), neighbors.end(),
            [&](NetworKit::node q) {
                return tree[q] == tree[p] && tree_capacity(q, p) > 0 && origin(q);
            });
        if (new_parent != neighbors.end()) {
            parent[p] = *new_parent;
            continue;
        }
        graph.forNeighborsOf(p, [&](NetworKit::node q) {
            if (tree[q] != tree[p]) {
                return;
            }
            if (tree_capacity(q, p) > 0) {
                active.push(q);
            }
            if (parent[q] == p) {
                parent[q] = NetworKit::none;
                orphan.push(q);
            }
        });
        tree[p] = NodeType::FREE;
    }
}

void BoykovKolmogorovMaximumFlow::run() {
    GraphTools::ensureDirectedGraph(graph);
    std::queue<NetworKit::node> active;
    std::queue<NetworKit::node> orphan;
    initialize(active);
    while (true) {
        const auto middle = grow(active);
        if (!middle) {
            break;
        }
        flow_size += augment(*middle, orphan);
        adopt(active, orphan);
    }
    hasRun = true;
}
}  // namespace Koala
