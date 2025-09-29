#include <flow/minimum_cost_flow/OrlinMCF.hpp>
#include <flow/KingRaoTarjanMaximumFlow.hpp>
#include <shortest_path/Dijkstra.hpp>

namespace Koala {

using edgeid = NetworKit::edgeid;
using node = NetworKit::node;

void OrlinMCF::init_flow() {
    this->graph.indexEdges();
    make_connected();
    make_uncapacitated();
    make_costs_nonnegative();
}

void OrlinMCF::make_costs_nonnegative() {
    // TODO 
}

void OrlinMCF::init_circulation() {
    construct_flow_from_circulation();
    init_flow();
}

void OrlinMCF::initialize() {
    flow.clear();
    original_graph = graph;
    excess = b;
    for(auto [key, value] : excess) {
        delta = std::max(delta, value); 
    }
}

long long OrlinMCF::cp(node u, node v) {
    return flow[{u, v}] - potential[u] + potential[v];
}

bool OrlinMCF::is_imbalanced() {
    for (auto [node, value] : excess) {
        if (value != 0) return true;
    }
    return false;
}

NetworKit::Graph OrlinMCF::generate_graph_for_sp(std::map<edge, std::pair<int,node>>& shortcuts) {
    NetworKit::Graph distGraph(graph.upperNodeIdBound(), true, true);
    auto [first, last] = uncapacitated_nodes_bounds;

    
    for (node i = first; i <= last; ++i) {
        if (graph.hasNode(i)) {
            node neighs[2];                
            for (int j = 0; j < 2; ++j) {
                neighs[j] = graph.getIthInNeighbor(i, j);
            }

            for (int j = 0; j < 2; ++j) {
                long long f = flow[{neighs[j], i}];

                if (f > 0) {
                    edge e = {neighs[j^1], neighs[j]};
                    long long weight = cp(e.first, i);
                    if (shortcuts.find(e) == shortcuts.end()) {
                        shortcuts[e] = {weight, i};
                    } else {
                        shortcuts[e] = std::min(shortcuts[e], {weight, i});
                    }

                    distGraph.addEdge(
                        e.first,
                        e.second,
                        weight
                    );
                }
            }
        }
    }
    graph.forEdges([&](node u, node v) {
        if (first <= v && v <= last) return;
        long long cost = cp(u, v);
        distGraph.addEdge(u, v, cost);
        if (flow[{u, v}] > 0) {
            distGraph.addEdge(v, u, -cost);
        }
    });
    return distGraph;
}

void OrlinMCF::fill_distances_uncapacitated(std::vector<double>& distances) {
    auto [first, last] = uncapacitated_nodes_bounds;   
    for(node i = first; i <= last; ++i) {
        if (!graph.hasNode(i)) 
            continue;

        graph.forInEdgesOf(i, [&](node _, node j) {
            distances[i] = std::min(distances[i], distances[j] + costs[{j, i}]);
        });
    }
}

void OrlinMCF::scaling_augmentation(node u, node v) {
    std::map<edge, std::pair<int, node>> shortcuts;
    NetworKit::Graph distGraph = generate_graph_for_sp(shortcuts);
    
    Dijkstra<FibonacciHeap> shortestPath(distGraph, u, true);
    shortestPath.run();
    std::vector<double> distances = shortestPath.getDistances();
    fill_distances_uncapacitated(distances);
    
    for(auto & [nd, pt] : potential) {
        pt -= lround(distances[nd]);
    }

    std::vector<node> path;

    if (is_added_uncapacitated(v)) {
        std::pair<long long, node> pathEnd = {INT64_MAX, -1};
        graph.forInEdgesOf(v, [&](node _, node w) {
            pathEnd = std::min(pathEnd, {distances[w] + costs[{w, v}], w}); 
        });
        node ending = pathEnd.second;
        path = shortestPath.getPath(ending);
        path.push_back(v);
    } else {
        path = shortestPath.getPath(v);
    }

    for (int i=0; i<path.size()-1; i++) {
        node p = path[i], q = path[i+1];

        if(flow[{q, p}] > 0) {
            flow[{q, p}] -= delta;
        } else {
            int e1, shortcut;
            e1 = shortcut = INT32_MAX;
            if (graph.hasEdge(p, q)) {
                e1 = costs[{p, q}];
            }
            auto iter = shortcuts.find({p, q});
            if (iter != shortcuts.end()) {
                shortcut = iter->second.first;
            }
            
            if (shortcut > e1) { 
                flow[{p, q}] += delta;
            } else {
                node x = iter->second.second;
                flow[{p, x}] += delta;
                flow[{q, x}] -= delta;
            }
        }
    }
    excess[path.front()] -= delta;
    excess[path.back()] += delta;        
}

void OrlinMCF::run_impl() {
    initialize();

    while (is_imbalanced()) {
        bool zero = true; 
        for (auto f : flow) {
            if (f.second) {
                zero = false;
                break;
            } 
        }
        long long newDelta{1};
        if (zero) { 
            long long maxDelta{0};
            for (auto ex : excess) {
                maxDelta = std::max(ex.second, maxDelta);
            }
            while (newDelta < maxDelta) newDelta <<= 1;
        }
        delta = std::min(delta, newDelta);

        delta_scaling_phase();
        delta /= 2;
    }

    uncontract_and_compute_flow();
}

void OrlinMCF::uncontract_and_compute_flow() {
    while (contractions.size()) {
        auto [u, v] = contractions.top();
        contractions.pop(); 
        potential[v] = potential[u];
    }

    NetworKit::Graph flowGraph(original_graph.numberOfNodes(), true, true);
    node s = flowGraph.addNode();
    node t = flowGraph.addNode();
    long long bsum{0};
    original_graph.forNodes([&](node n) { 
        long long bval = b[n];
        if (bval > 0) {
            flowGraph.addEdge(s, n, bval);
            bsum += bval;
        } else if (bval < 0) {
            flowGraph.addEdge(n, t, -bval);
        }
    });

    original_graph.forEdges([&](node u, node v){
        if (cp(u, v) == 0) {
            flowGraph.addEdge(u, v, bsum);
        }
    });

    KingRaoTarjanMaximumFlow maxflow(flowGraph, s, t);
    maxflow.run();
    original_graph.forEdges([&](node u, node v){
        flow[{u, v}] = maxflow.get_flow({u, v});
    });
}


void OrlinMCF::delta_scaling_phase() {
    contraction_phase();

    std::vector<NetworKit::node> s,t;
    graph.forNodes([&](NetworKit::node node) {
        if (excess[node] >= ALPHA*delta) {
            s.push_back(node);
        }
        if (excess[node] <= -ALPHA*delta) {
            t.push_back(node);
        }
    });

    while (s.size() && t.size()) {
        NetworKit::node k = s.back();
        NetworKit::node v = t.back();

        scaling_augmentation(k,v);
    }
}

void OrlinMCF::contraction_phase() {
    auto contractible = [&](edge& edge) -> bool {
        for (auto [e, f] : flow) {
            if (f >= 3*graph.numberOfNodes()*delta) {
                edge = e;
                return true;
            }
        }
        return false;
    };
    edge e;
    while(contractible(e)) {
        auto [u, v] = e;
        contract_nodes(u, v);
    }
}

void OrlinMCF::contract_nodes(node u, node v) {
    if (v < u) {
        std::swap(u, v);
    }

    std::set<edge> edge_set;
    graph.forEdgesOf(u, [&](node i, node j) {
        edge_set.insert({i, j});
    });
    graph.forInEdgesOf(u, [&](node i, node j) {
        edge_set.insert({j, i});
    });

    graph.forEdgesOf(v, [&](node, node w) {
        if(u != w) {
            if (edge_set.find({u, w}) == edge_set.end()) {
                edge_set.insert({u, w});
                graph.addEdge(u, w);
                flow[{u, w}] = flow[{v, w}];
                costs[{u, w}] = costs[{v, w}];
            } else {
                flow[{u, w}] += flow[{v, w}];
                costs[{u, w}] = std::min(costs[{u, w}], costs[{v, w}]);
            }
        }
        flow.erase({v, w});
        costs.erase({v, w});
    });
    graph.forInEdgesOf(v, [&](node, node w) {
        if (w != u) {
            if (edge_set.find({w, u}) == edge_set.end()) {
                edge_set.insert({w, u});
                graph.addEdge(w, u);
                flow[{w, u}] = flow[{w, v}];
                costs[{w, u}] = costs[{w, v}];
            } else {
                flow[{w, u}] += flow[{w, v}];
                costs[{w, u}] = std::min(costs[{w, u}], costs[{w, v}]);
            }
        }             
        flow.erase({w, v});
        costs.erase({w, v});
    });
    
    excess[u] = excess[u] + excess[v];
    excess.erase(v);

    graph.removeNode(v);
    contractions.push({u,v});
}

} /* namespace Koala */