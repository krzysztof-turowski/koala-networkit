#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>
#include <climits>
namespace Koala {

using edgeid = NetworKit::edgeid;
using node = NetworKit::node;
using Graph = NetworKit::Graph;

static void print_flows(edge_map<int> &printable);
static void forEachArcOf(Graph const&, node u, std::function<void(node, node, edgeid, bool)> f);

SuccessiveApproxMCC::SuccessiveApproxMCC(NetworKit::Graph& g,
    edge_map<MCFEdgeParams>& ep, node_map<int>& excess)
    : MinimumCostFlow(g,ep,excess) {
    constructCirculation();
}

SuccessiveApproxMCC::SuccessiveApproxMCC(NetworKit::Graph& g, 
    edge_map<MCFEdgeParams>& ep, NetworKit::node s, NetworKit::node t, int flow)
    : MinimumCostFlow(g, ep, s, t, flow) {
    constructCirculation();
}

inline double SuccessiveApproxMCC::cp(node v,
    node w, edgeid eid) {
    return static_cast<double>(costs[eid]) - potential[v] + potential[w];
}

inline int SuccessiveApproxMCC::uf(edgeid eid, bool reversed = false) {
    if(!reversed)
        return upperbound[eid] - flow[eid];
    return flow[eid] - lowerbound[eid];
}

void SuccessiveApproxMCC::force_flow(node v, node w,
    edgeid eid, int f) {
    flow[eid] += f;
    excess[v] -= f;
    excess[w] += f;
}

void SuccessiveApproxMCC::push(node v, node w, edgeid eid, bool reverse = false) {
    if (!reverse) { 
        int residualCapacity = upperbound[eid] - flow[eid];
        int change = std::min(residualCapacity, excess[v]);
        force_flow(v, w, eid, change);
    } else {
        int residualCapacity = flow[eid] - lowerbound[eid];
        int change = std::min(residualCapacity, excess[w]);
        force_flow(v, w, eid, -change);
    }
}

void SuccessiveApproxMCC::relabel(NetworKit::node const& v) {
    double mi = std::numeric_limits<double>::max();
    graph.forEdgesOf(
        v,
        [&](node i, node j, edgeid eid) {
            if (uf(eid) > 0) {
                mi = std::min(mi, potential[j]
                    + epsi
                    + costs[eid]);
            }
        });
    graph.forInEdgesOf(
        v,
        [&](node i, node j, edgeid eid) {
            if (uf(eid, true) > 0) {
                mi = std::min(mi, potential[j]
                    + epsi
                    - costs[eid]);
            }
        });
    potential[v] = mi;
}

void SuccessiveApproxMCC::refine() {
    epsi /= 2;
    graph.forEdges([&](node v, node w, edgeid eid) {
        double reducedCost = cp(v, w, eid);
        if (reducedCost < 0) {
            int residualCap = upperbound[eid] - flow[eid];
            force_flow(v, w, eid, residualCap);
        } else if (reducedCost > 0) {
            int residualCap = flow[eid] - lowerbound[eid];
            force_flow(v, w, eid, -residualCap); 
        }
    });

    wave();
}

void SuccessiveApproxMCC::wave() {
    DischargeList* list = new ToposortList(*this);

    NetworKit::node v = list->getNext();
    while (active.size()) {
        if (excess[v] > 0) {
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
    int& ex = excess[u];
    graph.forEdgesOf(u, [&](node u, node v, edgeid eid) {
        if (ex && cp(u, v, eid) < 0 && uf(eid) > 0) {
            push(u, v, eid);
        }
    });
    if (ex > 0) {
        // reversed arcs, that enters u can be outgoing in residual graph
        graph.forInEdgesOf(u, [&](node u, node v, node eid) {
            if (ex && cp(u, v, eid) > 0 && uf(eid, true) > 0) {
                push(u, v, eid, true);
            }
        });
    }
    if (ex > 0) {
        relabel(u);
        return true;
    }

    return false;
}

void SuccessiveApproxMCC::initialize() {
    potential.clear();
    excess.clear();
    flow.clear();
    active.clear();
    pr_id.clear();
    feasible = true;
    // std::cerr<<"INIT"<<"\n";
    for (auto [v,e] : excess) {
        if (e > 0) {
            active.insert(v);
        }
    }

    int mx = 0;

    for (auto [_, cost] : costs) {
        mx = std::max(mx, abs(cost));
    }
    epsi = static_cast<double>(mx);

    graph.forEdges([&](node u, node v, edgeid eid) {
        int cap;
        if ((cap = upperbound[eid]) < 0) {
            force_flow(u, v, eid, cap);
        } else if ((cap = lowerbound[eid] > 0)) {
            force_flow(u, v, eid, cap);
        }
    });
}

void SuccessiveApproxMCC::runImpl() {
    initialize();

    while (epsi >= 1.0/graph.numberOfNodes()) {
        refine();
        print_flows(flow);
    }

    min_cost = 0;
    
    graph.forEdges([&](node u, node v, edgeid eid) {
        min_cost += flow[eid] * costs[eid];
    });
}


SuccessiveApproxMCC::ToposortList::ToposortList(
    SuccessiveApproxMCC &approx) : approx(approx) {
    vis.clear();
    for (auto v : approx.graph.nodeRange()) {
        if(!vis[v]) dfs(v);
    }
    it2 = nodes.begin();
}

void SuccessiveApproxMCC::ToposortList::dfs(NetworKit::node v){
    vis[v] = true;

    forEachArcOf(approx.graph, v, [&](node i, node j, edgeid id, bool reversed) {
        if (approx.cp(i, j, id)*(reversed ? -1 : 1)  < 0 && approx.uf(id, reversed)) {
            if (!vis[j]) {
                dfs(j);
            }
        }
    });
    
    nodes.push_front(v);
}

NetworKit::node SuccessiveApproxMCC::ToposortList::getNext() {
    // std::cerr<<"NEXT\n";
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
        
        NetworKit::node value = *it1;
        nodes.erase(it1);
        nodes.push_front(value);
    }
}

static void print_flows(edgeid_map<int> &printable) {
    // std::cerr<< " --- AFTER REFINE TEST --- \n"; 
    for(auto [e, f] : printable){
        // std::cerr << "ARC {" << e.first<< ", "<<e.second<<"} FLOW: " << f <<  "\n";   
    }
    // std:: cerr<< " --- PRINT END --- \n\n";
}

static void forEachArcOf(
    Graph const& g,
    node u,
    std::function<void(node, node, edgeid, bool)> f) {
    
    g.forEdgesOf(u, [&](node u, node v, edgeid id) {
        f(u, v, id, false);
    });

    g.forInEdgesOf(u, [&](node u, node v, edgeid id) {
        f(u, v, id, true);
    });
}


} /* namespace Koala */
