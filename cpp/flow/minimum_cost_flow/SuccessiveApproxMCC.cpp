#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>
#include <climits>
// #define DEBUG_DUMP
#ifdef DEBUG_DUMP
#define DBG(x) std::cerr << x
#else 
#define DBG(x)
#endif

namespace Koala {

using node = NetworKit::node;

inline double SuccessiveApproxMCC::cp(uint64_t eid) {
    const Edge& edge = edges[eid];
    return static_cast<double>(edge.cost) - potential[edge.from] + potential[edge.to];
}

inline int64_t SuccessiveApproxMCC::uf(uint64_t eid) {
    return edges[eid].capacity - edges[eid].flow;
}

bool SuccessiveApproxMCC::is_imbalanced() {
    for (auto e : excess) {
        if (e) return true;
    }
    return false;
}

void SuccessiveApproxMCC::force_flow(uint64_t eid, int64_t f) {
    Edge& edge = edges[eid];
    edge.flow += f;
    edges[eid ^ 1].flow -= f;
    excess[edge.from] -= f;
    excess[edge.to] += f;
}

void SuccessiveApproxMCC::push(uint64_t eid) {
    node u = edges[eid].from;
    if (excess[u] > 0) {
        force_flow(eid, std::min(uf(eid), excess[u]));
    }
}

void SuccessiveApproxMCC::relabel(NetworKit::node const& u) {
    double mi = std::numeric_limits<double>::max();

    for (uint64_t eid : neigh_list[u]) {
        if (uf(eid) > 0) {
            const Edge& edge = edges[eid];
            mi = std::min(mi, potential[edge.to] + epsi + edge.cost);
        }
    }

    potential[u] = mi;
}

void SuccessiveApproxMCC::refine() {
    epsi /= 2;
    for (uint64_t eid = 0; eid < edges.size(); ++eid) {
        double reduced = cp(eid);
        if (reduced < 0) {
            force_flow(eid, uf(eid));
        }
    }
    wave();
}



void SuccessiveApproxMCC::wave() {
    DischargeList* list = new ToposortList(*this);

    NetworKit::node v = list->getNext();
    while (is_imbalanced()) {
        if (excess[v] > 0) {
            DBG("Discharging " << v << " with excess " << excess[v] << "\n");
            bool relabeled = discharge(v);
            if (relabeled) {
                list->moveToStart();
            }
        }
        v = list->getNext();
    }
    delete list;
}

bool SuccessiveApproxMCC::discharge(NetworKit::node const& u) {
    int64_t& ex = excess[u];
    for (uint64_t eid : neigh_list[u]) {
        DBG("Trying to push on edge " << edges[eid].from << "->" << edges[eid].to << " with cost " << edges[eid].cost << " and flow " << edges[eid].flow << "/" << edges[eid].capacity << "\n");
        DBG("Reduced cost: " << cp(eid) << " Residual capacity: " << uf(eid) << "\n");
        if (ex && cp(eid) < 0 && uf(eid) > 0) {
            push(eid);
        }
    }

    if (ex > 0) {
        DBG("Relabeling " << u << "\n");
        relabel(u);
        return true;
    }

    return false;
}

void SuccessiveApproxMCC::initialize() {
    auto& graph = network.getGraph();
    uint32_t nodeBound = graph.upperNodeIdBound();
    nodes_number = graph.numberOfNodes();
    potential.clear();
    excess.assign(nodeBound, 0);
   
    for (auto [key, value] : network.excess) {
        excess[key] = value;
    }

    potential.assign(nodeBound, 0);
    edges.reserve(2 * graph.numberOfEdges());
    int ptr = 0;
    
    neigh_list.assign(nodeBound, std::vector<uint64_t>());
    std::set<std::pair<NetworKit::node, NetworKit::node>> present;
    int64_t maxCost = 0;


    graph.forNodes([&](node u) { 
        graph.forNeighborsOf(u, [&](node v) {
            node x = std::min(u, v);
            node y = std::max(u, v);
            if (present.find({x, y}) != present.end())
                return;
            present.insert({x, y});

            uint32_t from = static_cast<uint32_t>(u);
            uint32_t to = static_cast<uint32_t>(v);
            int64_t cost = network.cost[{u, v}] - network.cost[{v, u}];
            int64_t capacity = network.capacity[{u, v}];
            int64_t capacity2 = network.capacity[{v, u}];

            neigh_list[from].push_back(edges.size());

            int64_t fl = std::min(capacity, (int64_t)0);
            fl = capacity2 < 0 ? -capacity2 : fl;

            edges.push_back({
                from, to,
                cost, capacity, fl
            });

            maxCost = std::max(maxCost, std::abs(cost));

            neigh_list[to].push_back(edges.size());
            edges.push_back({
                to, from,
                -cost, capacity2, -fl
            });

        });
    });

    epsi = static_cast<double>(maxCost);
}

void SuccessiveApproxMCC::run_impl() {
    initialize();

    while (epsi >= 1.0/nodes_number) {
        DBG("Refine with epsilon = " << epsi << "\n");
        refine();
    }

    min_cost = 0;
    
    for (const Edge& edge : edges) {
        min_cost += edge.flow * edge.cost;
        computed_flow[{edge.from, edge.to}] = edge.flow;
    }
    min_cost /= 2;
}


SuccessiveApproxMCC::ToposortList::ToposortList(
    SuccessiveApproxMCC &approx) : approx(approx) {
    vis.assign(approx.nodes_number, 0);
    auto& graph = approx.network.getGraph();
    for (auto v : graph.nodeRange()) {
        if(!vis[v]) dfs(v);
    }
    it2 = nodes.begin();
}

void SuccessiveApproxMCC::ToposortList::dfs(NetworKit::node u){
    vis[u] = true;

    for (auto eid : approx.neigh_list[u]) {
        if (approx.cp(eid) < 0 && approx.uf(eid) > 0) {
            auto v = approx.edges[eid].to;
            if (!vis[v]) {
                dfs(v);
            }
            
        }
    }

    nodes.push_front(u);
}

NetworKit::node SuccessiveApproxMCC::ToposortList::getNext() {
    if (it2 != nodes.end()) {
        it1 = it2;
        it2++;
    } else {
        it2 = nodes.begin();
        it1 = it2++;
    }
    return *it1;
}

void SuccessiveApproxMCC::ToposortList::moveToStart() {
    if (it1 != nodes.end()) {
        nodes.splice(nodes.begin(), nodes, it1);
    }
}

int64_t SuccessiveApproxMCC::getFlow(const NetworKit::Edge& edge) {
    return computed_flow[edge];
}

// static void print_flows(edgeid_map<int> &printable) {
//     // std::cerr<< " --- AFTER REFINE TEST --- \n"; 
//     for(auto [e, f] : printable){
//         // std::cerr << "ARC {" << e.first<< ", "<<e.second<<"} FLOW: " << f <<  "\n";   
//     }
//     // std:: cerr<< " --- PRINT END --- \n\n";
// }

} /* namespace Koala */
