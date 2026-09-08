#include <flow/minimum_cost_flow/SuccessiveApproximationMinimumCostFlow.hpp>

#include <algorithm>
#include <climits>
#include <limits>
#include <memory>
#include <vector>

namespace Koala {

using node = NetworKit::node;

inline double SuccessiveApproximationMinimumCostFlow::cp(NetworKit::index eid) {
    const Edge& edge = edges[eid];
    return static_cast<double>(edge.cost) - potential[edge.from] + potential[edge.to];
}

inline int64_t SuccessiveApproximationMinimumCostFlow::uf(NetworKit::index eid) {
    return edges[eid].capacity - edges[eid].flow;
}

bool SuccessiveApproximationMinimumCostFlow::is_imbalanced() {
    for (auto e : excess) {
        if (e) return true;
    }
    return false;
}

void SuccessiveApproximationMinimumCostFlow::force_flow(NetworKit::index eid, int64_t f) {
    Edge& edge = edges[eid];
    edge.flow += f;
    edges[eid ^ 1].flow -= f;
    excess[edge.from] -= f;
    excess[edge.to] += f;
}

void SuccessiveApproximationMinimumCostFlow::push(NetworKit::index eid) {
    node u = edges[eid].from;
    if (excess[u] > 0) {
        force_flow(eid, std::min(uf(eid), excess[u]));
    }
}

void SuccessiveApproximationMinimumCostFlow::relabel(NetworKit::node const& u) {
    double mi = std::numeric_limits<double>::infinity();

    for (NetworKit::index eid : neighbors[u]) {
        if (uf(eid) > 0) {
            const Edge& edge = edges[eid];
            mi = std::min(mi, potential[edge.to] + epsi + edge.cost);
        }
    }

    potential[u] = mi;
}

void SuccessiveApproximationMinimumCostFlow::refine() {
    epsi /= 2;
    for (NetworKit::index eid = 0; eid < edges.size(); ++eid) {
        double reduced = cp(eid);
        if (reduced < 0) {
            force_flow(eid, uf(eid));
        }
    }
    wave();
}

void SuccessiveApproximationMinimumCostFlow::wave() {
    std::unique_ptr<DischargeList> list = std::make_unique<ToposortList>(*this);

    NetworKit::node v = list->getNext();
    while (is_imbalanced()) {
        if (excess[v] > 0) {
            bool relabeled = discharge(v);
            if (relabeled) {
                list->moveToStart();
            }
        }
        v = list->getNext();
    }
}

bool SuccessiveApproximationMinimumCostFlow::discharge(NetworKit::node const& u) {
    int64_t& ex = excess[u];
    for (NetworKit::index eid : neighbors[u]) {
        if (ex && cp(eid) < 0 && uf(eid) > 0) {
            push(eid);
        }
    }

    if (ex > 0) {
        relabel(u);
        return true;
    }

    return false;
}

void SuccessiveApproximationMinimumCostFlow::initialize() {
    auto& graph = network.getGraph();
    NetworKit::count nodeBound = graph.upperNodeIdBound();
    nodes_number = graph.numberOfNodes();
    potential.clear();
    excess.assign(nodeBound, 0);

    for (auto [key, value] : network.excess) {
        excess[key] = value;
    }

    potential.assign(nodeBound, 0);
    edges.reserve(2 * graph.numberOfEdges());

    neighbors.assign(nodeBound, std::vector<NetworKit::index>());
    int64_t max_cost = 0;

    graph.forNodes([&](node u) {
        graph.forNeighborsOf(u, [&](node v) {
            node from = u;
            node to = v;
            int64_t cost = network.cost[{u, v}];
            int64_t capacity = network.capacity[{u, v}];

            neighbors[from].push_back(edges.size());
            edges.push_back({
                from, to,
                cost, capacity, 0LL
            });

            max_cost = std::max(max_cost, std::abs(cost));

            neighbors[to].push_back(edges.size());
            edges.push_back({
                to, from,
                -cost, 0LL, 0LL
            });
        });
    });

    epsi = static_cast<double>(max_cost);
}

void SuccessiveApproximationMinimumCostFlow::run_impl() {
    initialize();

    while (epsi >= 1.0/nodes_number) {
        refine();
    }

    min_cost = 0;

    for (const Edge& edge : edges) {
        min_cost += edge.flow * edge.cost;
        computed_flow[{edge.from, edge.to}] = edge.flow;
    }
    min_cost /= 2;
}

SuccessiveApproximationMinimumCostFlow::ToposortList::ToposortList(
    SuccessiveApproximationMinimumCostFlow &approx) : approx(approx) {
    vis.assign(approx.nodes_number, 0);
    auto& graph = approx.network.getGraph();
    for (auto v : graph.nodeRange()) {
        if (!vis[v]) dfs(v);
    }
    it2 = nodes.begin();
}

void SuccessiveApproximationMinimumCostFlow::ToposortList::dfs(NetworKit::node u) {
    vis[u] = true;

    for (auto eid : approx.neighbors[u]) {
        if (approx.cp(eid) < 0 && approx.uf(eid) > 0) {
            auto v = approx.edges[eid].to;
            if (!vis[v]) {
                dfs(v);
            }
        }
    }

    nodes.push_front(u);
}

NetworKit::node SuccessiveApproximationMinimumCostFlow::ToposortList::getNext() {
    if (it2 != nodes.end()) {
        it1 = it2;
        it2++;
    } else {
        it2 = nodes.begin();
        it1 = it2++;
    }
    return *it1;
}

void SuccessiveApproximationMinimumCostFlow::ToposortList::moveToStart() {
    if (it1 != nodes.end()) {
        nodes.splice(nodes.begin(), nodes, it1);
    }
}

int64_t SuccessiveApproximationMinimumCostFlow::getFlow(const NetworKit::Edge& edge) {
    return computed_flow[edge];
}


} /* namespace Koala */
