#include <flow/minimum_cost_flow/SuccessiveApproxMCF.hpp>
#include <climits>
namespace Koala {

inline int SuccessiveApproxMCF::cp(NetworKit::node const& v, NetworKit::node const& w){
    return edge_params[{v,w}].cost - price[v] + price[w];
} 

inline int SuccessiveApproxMCF::uf(NetworKit::node const& v, NetworKit::node const& w){
    return edge_params[{v,w}].capacity - flow[{v,w}];
}

void SuccessiveApproxMCF::force_flow(NetworKit::node const& v, NetworKit::node const& w, int f){
    flow[{v,w}]+=f;
    flow[{w,v}]-=f;
    excess[v] -= f;
    excess[w] += f;

    if(excess[v] > 0) active.insert(v);
    else active.erase(v);
    if(excess[w] > 0) active.insert(w);
    else active.erase(w);
}



void SuccessiveApproxMCF::push(NetworKit::node const& v, NetworKit::node const& w){
    MCFEdgeParams edge = edge_params[{v, w}];
    int res = edge.capacity - flow[{v,w}];
    int change = std::min(res, excess[v]);
    force_flow(v,w, change);
}

void SuccessiveApproxMCF::relabel(NetworKit::node const& v){
    double mi = INFINITY;
    graph.forEdgesOf(
        v,
        [&](NetworKit::node V, NetworKit::node W){
            if(uf(V,W) > 0){
                mi = std::min(mi,price[W] 
                    + edge_params[{V,W}].capacity 
                    + epsi);
            }
        }
    );
    price[v] = mi;
}

void SuccessiveApproxMCF::push_relabel(NetworKit::node const& v){
    int id = pr_id[v];
    NetworKit::node w = graph.getIthNeighbor(v, id);
    if(uf(v,w) > 0 && cp(v,w) < 0){
        push(v,w);
    }
    else{
        if(graph.degreeOut(v) != id+1){
            pr_id[v]++;
        }
        else{
            pr_id[v]=0;
            relabel(v);
        }
    }
}

void SuccessiveApproxMCF::refine(){
    epsi /= 2;
    graph.forEdges(
        [&](NetworKit::node v, NetworKit::node w) {
            if(cp(v,w) < 0){
                int add = edge_params[{v,w}].capacity - flow[{v,w}];
                force_flow(v,w,add);
            }
        }
    );

    while(active.size()){
        push_relabel(*active.begin());
    }
}

void SuccessiveApproxMCF::initialize(){
    price.clear();
    flow.clear();
    active.clear();
    excess.clear();
    pr_id.clear();
    feasible = true;

    int mx = 0;
    for(auto [a,b] : edge_params){
        if(int c = b.capacity < 0){
            force_flow(a.first, a.second, c);
        }
        mx = std::max(mx, abs(b.cost));
    }
    epsi = (double)mx;
}

void SuccessiveApproxMCF::runImpl(){
    initialize();

    if(epsi >= 1.0/graph.numberOfNodes()){
        refine();
        for(auto [v, e] : excess){
            if(e > 0){
                feasible = false;
                return;
            }        
        }
    }


    while(epsi >= 1.0/graph.numberOfNodes()){
        refine();
    }

    min_cost = 0;
    graph.forEdges(
        [&](NetworKit::node v, NetworKit::node w) {
            min_cost += flow[{v,w}]*edge_params[{v,w}].cost;
        }
    );
    min_cost/=2;
}

}