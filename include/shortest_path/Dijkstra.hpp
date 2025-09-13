/*
 * Dijkstra.hpp
 *  
 * Implementation borrowed from the NetworKit's Dijkstra algorithm 
 * and extended to allow using own heaps. 
 * 
 * 
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 *  Modified on: Aug 29, 2025
 *      Author: Oskar Krygier 
 */

#pragma once

#include <networkit/distance/SSSP.hpp>

#include <structures/heap/FibonacciHeap.hpp>

namespace Koala {
    
template<template <typename...> class THeap> 
class Dijkstra final : public NetworKit::SSSP {
 private:
    THeap<std::pair<double, NetworKit::node>, std::greater<std::pair<double, NetworKit::node>>> heap;

 public:
    Dijkstra(const NetworKit::Graph &G, NetworKit::node source,
                    bool storePaths = false, bool storeNodesSortedByDistance = false,
                    NetworKit::node target = NetworKit::none)
        : SSSP(G, source, storePaths, storeNodesSortedByDistance, target) {}

    void run() override;
};

template <template <typename...> class THeap> 
void Dijkstra<THeap>::run() {
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
    heap.push({0., source});

    auto initPath = [&](NetworKit::node u, NetworKit::node v) {
        if (storePaths) {
            previous[v] = {u};
            npaths[v] = npaths[u];
        }
    };
    bool breakWhenFound = (target != NetworKit::none);

    do {
        // TRACE("pq size: ", heap.size());
        auto [_, u] = heap.top();
        // std::cerr<<"new guy " << u << "\n";
        heap.pop();
        
        sumDist += distances[u];

        if ((breakWhenFound && target == u) || distances[u] == infDist)
            break;

        if (storeNodesSortedByDistance)
            nodesSortedByDistance.push_back(u);

        G->forNeighborsOf(u, [&](NetworKit::node v, NetworKit::edgeweight w) {
            // std::cerr << u << " " << v << ": " << w <<'\n';
            double newDist = distances[u] + w;
            if (distances[v] == infDist) {
                distances[v] = newDist;
                heap.push({distances[v], v});
                ++reachedNodes;
                if (storePaths)
                    initPath(u, v);
            } else if (distances[v] > newDist) {
                if (storePaths)
                    initPath(u, v);
                distances[v] = newDist;
                heap.push({distances[v], v});
            } else if (storePaths && distances[v] == newDist) {
                previous[v].push_back(u);
                npaths[v] += npaths[u];
            }
        });
    } while (!heap.empty());

    hasRun = true;
}

} /* namespace Koala */