#pragma once

#include <networkit/distance/SSSP.hpp>

#include <structures/heap/FibonacciHeap.hpp>

namespace Koala {
    
class FibonacciDijkstra final : public NetworKit::SSSP {
 private:
    FibonacciHeap<NetworKit::node> heap;

 public:
    FibonacciDijkstra(const NetworKit::Graph &G, NetworKit::node source,
                    bool storePaths = false, bool storeNodesSortedByDistance = false,
                    NetworKit::node target = NetworKit::none)
        : SSSP(G, source, storePaths, storeNodesSortedByDistance, target) {}

    void run() override;
};

void FibonacciDijkstra::run() {
    auto infDist = std::numeric_limits<NetworKit::edgeweight>::max();
    std::fill(distances.begin(), distances.end(), infDist);

    if (distances.size() < G->upperNodeIdBound())
        distances.resize(G->upperNodeIdBound(), infDist);

    sumDist = 0.;
    reachedNodes = 1;

    if (storePaths) {
        previous.clear();
        previous.resize(G->upperNodeIdBound());
        npaths.clear();
        npaths.resize(G->upperNodeIdBound(), 0);
        npaths[source] = 1;
    }

    if (storeNodesSortedByDistance) {
        nodesSortedByDistance.clear();
        nodesSortedByDistance.reserve(G->upperNodeIdBound());
    }

    // priority queue with distance-node pairs
    distances[source] = 0.;
    heap.clear();
    heap.push(source);

    auto initPath = [&](NetworKit::node u, NetworKit::node v) {
        if (storePaths) {
            previous[v] = {u};
            npaths[v] = npaths[u];
        }
    };
    bool breakWhenFound = (target != NetworKit::none);

    do {
        // TRACE("pq size: ", heap.size());
        NetworKit::node u = heap.top();
        sumDist += distances[u];

        if ((breakWhenFound && target == u) || distances[u] == infDist)
            break;

        if (storeNodesSortedByDistance)
            nodesSortedByDistance.push_back(u);

        G->forNeighborsOf(u, [&](NetworKit::node v, NetworKit::edgeweight w) {
            double newDist = distances[u] + w;
            if (distances[v] == infDist) {
                distances[v] = newDist;
                heap.push(v);
                ++reachedNodes;
                if (storePaths)
                    initPath(u, v);
            } else if (distances[v] > newDist) {
                if (storePaths)
                    initPath(u, v);
                distances[v] = newDist;
                heap.push(v);
            } else if (storePaths && distances[v] == newDist) {
                previous[v].push_back(u);
                npaths[v] += npaths[u];
            }
        });
    } while (!heap.empty());

    hasRun = true;
}

} /* namespace Koala */