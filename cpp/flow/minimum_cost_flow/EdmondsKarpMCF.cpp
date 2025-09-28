#include <flow/minimum_cost_flow/EdmondsKarpMCF.hpp>
#include <shortest_path/Dijkstra.hpp>
#include <vector>

using node = NetworKit::node;
using edgeid = NetworKit::edgeid;
using edgeweight = NetworKit::edgeweight;

namespace Koala {

void EdmondsKarpMCF::initialize() {
    excess = b;
    flow.clear();
    potential.clear();
    NetworKit::edgeweight maxWeight = 0;
    makeConnected();
    
    for (auto edge : graph.edgeWeightRange()) {
        maxWeight = std::max(maxWeight, edge.weight);
    }
    delta = 1;
    while (delta < lround(maxWeight)) delta <<= 1;
}

long long EdmondsKarpMCF::cp(node u, node v) {
    return costs[{u, v}] - potential[u] + potential[v];
}

void EdmondsKarpMCF::send(node u, node v, long long value) {
    excess[u] -= value;
    excess[v] += value;
    flow[{u, v}] += value;
}

NetworKit::Graph EdmondsKarpMCF::getDeltaResidual() {
    NetworKit::Graph deltaResidual(graph.numberOfNodes(), true, true);
    graph.forEdges([&](node u, node v, edgeweight weight) {
        long long reducedCost = cp(u, v);
        int f = flow[{u, v}];

        if (f + delta <= lround(weight)) {
            deltaResidual.addEdge(u, v, reducedCost);
        }
        if (f >= delta) {
            deltaResidual.addEdge(v, u, -reducedCost);
        }
    });
    return deltaResidual;
}

void EdmondsKarpMCF::deltaScalingPhase() {

    graph.forEdges([&](node u, node v, edgeweight weight) {
        long long f = flow[{u, v}];

        if(f >= delta && cp(u, v) > 0) {
            send(u, v, -f);
        }
        else if (f + delta <= lround(weight) && cp(u, v) < 0){ 
            send(u, v, lround(weight) - f);
        }
    });

    std::vector<node> S, T;
    graph.forNodes([&](node v){
        int ex = excess[v];
        if (ex >= delta) S.push_back(v);
        else if (ex <= -delta) T.push_back(v);
    });

    while (!S.empty() && !T.empty()) {
        auto deltaResidual = getDeltaResidual();
        node k = S.back();
        node l = T.back();

        // shortest path
        auto dijkstra = Dijkstra<FibonacciHeap>(deltaResidual, k, true);
        dijkstra.run();
        std::vector<node> path = dijkstra.getPath(l);
        std::vector<double> distances = dijkstra.getDistances();
        
        // update potential
        for(node v : graph.nodeRange()) {
            potential[v] -= lround(distances[v]);
        }

        // augment
        for (int i=0; i<path.size()-1; i++) {
            node p = path[i], q = path[i+1];

            if (flow[{q, p}] >= delta) {
                send(q, p, -delta);
            } else {
                send(p, q, delta);
            }
        }

        // update S, T
        if (excess[k] < delta) S.pop_back();
        if (excess[l] > -delta) T.pop_back();
    }
}

void EdmondsKarpMCF::runImpl() {
    initialize();

    while (delta >= 1) {
        deltaScalingPhase();
        delta /= 2;
    }

    min_cost = 0;
    graph.forEdges([&](node u, node v, edgeid eid) {
        min_cost += getFlow({u, v}) * costs[{u, v}];
    });
}

} /* namespace Koala */