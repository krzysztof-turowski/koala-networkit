#include <flow/minimum_cost_flow/OrlinMCF.hpp>
#include <flow/KingRaoTarjanMaximumFlow.hpp>
#include <shortest_path/Dijkstra.hpp>

namespace Koala {

using edgeid = NetworKit::edgeid;
using node = NetworKit::node;

void OrlinMCF::initFlow() {
    this->graph.indexEdges();
    makeConnected();
    makeUncapacitated();
    makeCostsNonnegative();
}

void OrlinMCF::makeCostsNonnegative() {
    
}

void OrlinMCF::initCirculation() {
    constructFlowFromCirculation();
    initFlow();
}

void OrlinMCF::initialize() {
    flow.clear();
    originalGraph = graph;
    excess = b;
    for(auto [key, value] : excess) {
        delta = std::max(delta, value); 
    }
}

int OrlinMCF::cp(node i, node j, edgeid id) {
    return flow[id] - potential[i] + potential[j];
}

bool OrlinMCF::isImbalanced() {
    for (auto [node, value] : excess) {
        if (value != 0) return true;
    }
    return false;
}

NetworKit::Graph OrlinMCF::generateGraphForSP() {
    NetworKit::Graph distGraph(graph.upperNodeIdBound(), true, true);
    auto [first, last] = uncapacitatedNodesBounds;
    for (node i = first; i <= last; ++i) {
        if (graph.hasNode(i)) {
            std::vector<std::pair<node, edgeid>> neighs;                
            graph.forInEdgesOf(i, 
                [&](node u, node v, edgeid id) {
                neighs.push_back({v,id});
            });

            for (int j = 0; j < 2; ++j) {
                int f = flow[neighs[j].second];

                if (f > 0) {
                    distGraph.addEdge(
                        neighs[j ^ 1].first,
                        neighs[j].first,
                        cp(neighs[j ^ 1].first, i, neighs[j ^ 1].second)
                    );
                }
            }
        }
    }
    graph.forEdges([&](node u, node v, edgeid id) {
        if (first <= v && v <= last) return;
        int cost = cp(u, v, id);
        distGraph.addEdge(u, v, cost);
        if (flow[id] > 0) {
            distGraph.addEdge(v, u, -cost);
        }
    });
    return distGraph;
}

void OrlinMCF::scalingAugmentation(node u, node v) {
    NetworKit::Graph distGraph = generateGraphForSP();
    
    Dijkstra<FibonacciHeap> shortestPath(distGraph, u, true);
    shortestPath.run();
    auto path = shortestPath.getPath(v);
    std::vector<double> distances = shortestPath.getDistances();
    for(auto & [nd, pt] : potential) {
        pt -= lround(distances[nd]);
    }
    std::vector<std::pair<edgeid, bool>> pathIds;
    for (int i=0; i<path.size()-1; i++) {
        std::tuple<int, edgeid, bool> mini = {INT32_MAX, -1, false};
        node p = path[i], q = path[i+1];
        graph.forEdgesOf(p, [&](node u, node v, edgeid id) {
            if (v != q) return;
            mini = std::min(mini, {costs[id], id , false});
        });
        graph.forInEdgesOf(p, [&](node u, node v, edgeid id) {
            if (v != q || flow[id] < delta) return;
            mini = std::min(mini, {costs[id], id, true});
        });
        auto [_, id, reversed] = mini;
        pathIds.push_back({id, reversed});
    }
    augment(v, u, pathIds, delta);
}

void OrlinMCF::runImpl() {
    initialize();

    while (isImbalanced()) {
        bool zero = true; 
        for (auto f : flow) {
            if (f.second) {
                zero = false;
                break;
            } 
        }
        int maxDelta = 0;
        if (zero) { 
            for (auto ex : excess) {
                maxDelta = std::max(ex.second, maxDelta);
            }
        }
        delta = std::min(delta, maxDelta);

        deltaScalingPhase();
        delta /= 2;
    }

    uncontractAndComputeFlow();
}

void OrlinMCF::uncontractAndComputeFlow() {
    while (contractions.size()) {
        auto [u, v] = contractions.top();
        contractions.pop(); 
        potential[v] = potential[u];
    }

    NetworKit::Graph flowGraph(originalGraph.numberOfNodes(), true, true);
    node s = flowGraph.addNode();
    node t = flowGraph.addNode();
    int bsum{0};
    originalGraph.forNodes([&](node n) { 
        int bval = b[n];
        if (bval > 0) {
            flowGraph.addEdge(s, n, bval);
            bsum += bval;
        } else if (bval < 0) {
            flowGraph.addEdge(n, t, -bval);
        }
    });

    originalGraph.forEdges([&](node u, node v, edgeid id){
        if (cp(u, v, id) == 0) {
            flowGraph.addEdge(u, v, bsum);
        }
    });

    KingRaoTarjanMaximumFlow maxflow(flowGraph, s, t);
    maxflow.run();
    originalGraph.forEdges([&](node u, node v, edgeid id){
        flow[id] = maxflow.get_flow({u, v});
    });
}


void OrlinMCF::deltaScalingPhase() {
    contractionPhase();

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

        scalingAugmentation(k,v);
    }
}

void OrlinMCF::augment(node u, node v, edgeid edgeId, bool isReversed, int f) {
    augment(u, v, {{edgeId, isReversed}}, f);
} 

void OrlinMCF::augment(node u, node v, const std::vector<std::pair<edgeid, bool>> &edgeIds, int f) {
    for(auto [edgeId, reversed] : edgeIds) {
        flow[edgeId] += (!reversed ? f : -f);
    }

    excess[u] -= f;
    excess[v] += f;
}

void OrlinMCF::contractionPhase() {
    auto contractible = [&](edgeid& id) -> bool {
        for (auto [eid, f] : flow) {
            if (f >= 3*graph.numberOfNodes()*delta) {
                id = eid;
                return;
            }
        }
    };
    edgeid id;
    while(contractible(id)) {
        auto [u, v] = graph.edgeById(id);
        contractNodes(u, v);
    }
}

void OrlinMCF::contractNodes(NetworKit::node v, NetworKit::node w) {
    if (w < v) {
        std::swap(v, w);
    }

    graph.forEdgesOf(w, 
        [&](NetworKit::node, NetworKit::node y, NetworKit::edgeid id) {
            if(y != v) {
                graph.addEdge(v, y);
                int upper = graph.upperEdgeIdBound();
                flow[upper - 1] = flow[id];
                costs[upper - 1] = costs[id];
            }
            
            flow.erase(id);
            costs.erase(id);
    });
    graph.forInEdgesOf(w,
        [&](NetworKit::node, NetworKit::node y, NetworKit::edgeweight weight, NetworKit::edgeid id) {
            if (y != v) {
                graph.addEdge(y, v, weight);
                int upper = graph.upperEdgeIdBound();
                flow[upper - 1] = flow[id];
                costs[upper - 1] = costs[id];
            }             
            flow.erase(id);
            costs.erase(id);
    });
    
    excess[v] = excess[w] + excess[v];
    excess.erase(w);
    // b[v] = b[v] + b[w];
    // b.erase(w);
    graph.removeNode(w);
    contractions.push({v,w});
}

} /* namespace Koala */