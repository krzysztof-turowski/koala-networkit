#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>
#include <climits>
namespace Koala {

static void print_flows(edge_map<int> &printable);

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

inline double SuccessiveApproxMCC::cp(NetworKit::node const& v,
    NetworKit::node const& w) {
    return static_cast<double>(edge_params[{v, w}].cost) - price[v] + price[w];
}

inline int SuccessiveApproxMCC::uf(NetworKit::node const& v,
    NetworKit::node const& w) {
    return edge_params[{v, w}].capacity - flow[{v, w}];
}

void SuccessiveApproxMCC::force_flow(NetworKit::node const& v,
    NetworKit::node const& w, int f) {
    flow[{v, w}] += f;
    flow[{w, v}] -= f;
    excess[v] -= f;
    excess[w] += f;

    if (excess[v] > 0)
        active.insert(v);
    else
        active.erase(v);

    if (excess[w] > 0)
        active.insert(w);
    else
        active.erase(w);
}

void SuccessiveApproxMCC::push(NetworKit::node const& v,
    NetworKit::node const& w) {
    MCFEdgeParams edge = edge_params[{v, w}];
    int res = edge.capacity - flow[{v, w}];
    int change = std::min(res, excess[v]);
    force_flow(v, w, change);
}

void SuccessiveApproxMCC::relabel(NetworKit::node const& v) {
    double mi = INFINITY;
    graph.forEdgesOf(
        v,
        [&](NetworKit::node V, NetworKit::node W) {
            if (uf(V, W) > 0) {
                mi = std::min(mi, price[W]
                    + epsi
                    + edge_params[{V, W}].cost);
            }
        });
    price[v] = mi;
}

void SuccessiveApproxMCC::push_relabel(NetworKit::node const& v) {
    // std::cerr << "PUSH RELABEL NODE " << v << "\n"; 
    while (excess[v] > 0) {
        int id = pr_id[v];
        NetworKit::node w = graph.getIthNeighbor(v, id);
        if (uf(v, w) > 0 && cp(v, w) < 0) {
            // std::cerr<<"PUSH {" << v<< ", " << w << "} << \n";
            push(v, w);
        } else {
            if (graph.degreeOut(v) != id+1) {
                pr_id[v]++;
            } else {
                // std::cerr<<"RELABEL {" << v << "}\n";

                pr_id[v] = 0;
                relabel(v);
            }
        }
    }
}

void SuccessiveApproxMCC::refine() {
    epsi /= 2;
    for(auto v : graph.nodeRange()){
        // std::cerr<<"NODE "<<v<<"\n";
        for(auto w : graph.neighborRange(v)){
            // std::cerr << "FOR EDGE {"<<v<<", "<<w<<"} cp equals " << cp(v,w)<<"\n";
            if (cp(v, w) < 0) {
                int add = edge_params[{v, w}].capacity - flow[{v, w}];
                force_flow(v, w, add);
            }
        }
    
    }

    while (active.size()) {
        push_relabel(*active.begin());
    }
}

void SuccessiveApproxMCC::initialize() {
    price.clear();
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
    for (auto [a, b] : edge_params) {
        if (int c = b.capacity < 0) {
            force_flow(a.first, a.second, c);
        }
        mx = std::max(mx, abs(b.cost));
    }
    epsi = static_cast<double>(mx);
}

void SuccessiveApproxMCC::runImpl() {
    initialize();
    // std::cerr << "RUN IMPL" <<"\n";
    if (epsi >= 1.0/graph.numberOfNodes()) {
        refine();
        for (auto [v, e] : excess) {
            if (e > 0) {
                feasible = false;
                return;
            }
        }
    }
    print_flows(flow);


    while (epsi >= 1.0/graph.numberOfNodes()) {
        refine();
        print_flows(flow);
    }

    min_cost = 0;
    
    for (auto v : graph.nodeRange()) {
        for (auto w : graph.neighborRange(v)) {
            min_cost += flow[{v, w}]*edge_params[{v, w}].cost;
            // std::cerr << "Add cost for {" << v<< ", "<<w<<"} with cost "<< edge_params[{v, w}].cost << " for flow " << flow[{v, w}] << "\n";
        }
    }

    min_cost/=2;
}

static void print_flows(edge_map<int> &printable) {
    // std::cerr<< " --- AFTER REFINE TEST --- \n"; 
    for(auto [e, f] : printable){
        // std::cerr << "ARC {" << e.first<< ", "<<e.second<<"} FLOW: " << f <<  "\n";   
    }
    // std:: cerr<< " --- PRINT END --- \n\n";
}

} /* namespace Koala */
