#include <flow/minimum_cost_flow/EdmondsKarpMCF.hpp>
#include <shortest_path/Dijkstra.hpp>
#include <vector>

using node = NetworKit::node;
using edgeid = NetworKit::edgeid;

namespace Koala {

void EdmondsKarpMCF::initialize() {
    // we assume that the network is uncapacitated OK NOT REALLY

    flow.clear();
    potential.clear();
    NetworKit::edgeweight maxWeight = 0;

    for (auto edge : graph.edgeWeightRange()) {
        maxWeight = std::max(maxWeight, edge.weight);
    }
    delta = 1;
    while (delta < maxWeight) delta <<= 1;
}

int EdmondsKarpMCF::cp(node u, node v, edgeid eid) {
    return costs[eid] - potential[u] + potential[v];
}

void EdmondsKarpMCF::send(node u, node v, edgeid eid, int value) {
    excess[u] -= value;
    excess[v] += value;
    flow[eid] += value;
}

NetworKit::Graph EdmondsKarpMCF::getDeltaResidual() {
    NetworKit::Graph deltaResidual(graph.numberOfNodes(), true, true);
    graph.forEdges([&](node u, node v, edgeid id) {
        int reducedCost = cp(u, v, id);
        deltaResidual.addEdge(u, v, reducedCost);
        if (flow[id] >= delta) {
            deltaResidual.addEdge(v, u, -reducedCost);
        }
    });
    return deltaResidual;
}

void EdmondsKarpMCF::deltaScalingPhase() {
    graph.forEdges([&](node u, node v, edgeid id) {
        int resCap = flow[id];
        if(resCap >= delta && cp(u, v, id) > 0) {
            send(u, v, id, -resCap);
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
            std::pair<int, edgeid> mini = {INT32_MAX, -1};
            node p = path[i], q = path[i+1];
            graph.forEdgesOf(p, [&](node u, node v, edgeid id) {
                if (v != q) return;
                mini = std::min(mini, {costs[id], id});
            });
            graph.forInEdgesOf(p, [&](node u, node v, edgeid id) {
                if (v != q || flow[id] < delta) return;
                mini = std::min(mini, {costs[id], id});
            });
            send(p, q, mini.second, delta);
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
}

} /* namespace Koala */