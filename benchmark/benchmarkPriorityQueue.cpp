#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <string>
#include <type_traits>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include <benchmark/utils.hpp>
#include "structures/PriorityQueue.hpp"
#include "structures/heap/RankPairingHeap.hpp"
#include "structures/heap/SkewHeap.hpp"
#include "structures/heap/WeakHeap.hpp"

class Dijkstra {
 private:
    NetworKit::Graph& G;
    NetworKit::node source;
    std::vector<double> distances;
    std::vector<NetworKit::node> predecessors;
    std::vector<bool> visited;

 public:
    Dijkstra(NetworKit::Graph& G, NetworKit::node source) : G(G), source(source) {
        distances.resize(G.upperNodeIdBound(), std::numeric_limits<double>::infinity());
        predecessors.resize(G.upperNodeIdBound(), NetworKit::none);
        visited.resize(G.upperNodeIdBound(), false);
    }

    template <class PQType>
    void run() {
        static_assert(std::is_base_of<Koala::PriorityQueue<NetworKit::node>, PQType>::value,
                     "PQType must inherit from Koala::PriorityQueue");

        PQType pq;

        distances[source] = 0;
        pq.push(source);

        while (!pq.empty()) {
            NetworKit::node u = pq.top();
            pq.pop();

            if (visited[u]) continue;

            visited[u] = true;

            G.forNeighborsOf(u, [&](NetworKit::node v, double weight) {
                double newDist = distances[u] + weight;

                if (newDist < distances[v]) {
                    distances[v] = newDist;
                    predecessors[v] = u;
                    pq.push(v);
                }
            });
        }
    }

    std::vector<double> getDistances() const {
        return distances;
    }

    double getDistance(NetworKit::node t) const {
        return distances[t];
    }

    std::vector<NetworKit::node> getPath(NetworKit::node t) const {
        std::vector<NetworKit::node> path;
        if (distances[t] == std::numeric_limits<double>::infinity()) {
            return path;
        }

        for (NetworKit::node v = t; v != source; v = predecessors[v]) {
            path.push_back(v);
            if (v == predecessors[v]) break;
        }
        path.push_back(source);
        std::reverse(path.begin(), path.end());
        return path;
    }
};

enum class Algorithm : uint32_t { RANK_PAIRING_HEAP, SKEW_HEAP, WEAK_HEAP };

std::map<std::string, Algorithm> ALGORITHM = {
    { "rank_pairing_heap", Algorithm::RANK_PAIRING_HEAP },
    { "skew_heap", Algorithm::SKEW_HEAP },
    { "weak_heap", Algorithm::WEAK_HEAP }
};

void runDijkstra(
        NetworKit::Graph& G, NetworKit::node source, const std::string& algorithm_name,
        Algorithm algorithm) {
    auto start = std::chrono::high_resolution_clock::now();

    Dijkstra dijkstra(G, source);

    switch (algorithm) {
        case Algorithm::RANK_PAIRING_HEAP:
            dijkstra.run<Koala::RankPairingHeap<NetworKit::node>>();
            break;
        case Algorithm::SKEW_HEAP:
            dijkstra.run<Koala::SkewHeap<NetworKit::node>>();
            break;
        case Algorithm::WEAK_HEAP:
            dijkstra.run<Koala::WeakHeap<NetworKit::node>>();
            break;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Algorithm: Dijkstra with " << algorithm_name << std::endl;
    std::cout << "Time: " << duration.count() << "s" << std::endl;

    double totalDist = 0;
    int reachableNodes = 0;
    for (NetworKit::node v = 0; v < G.upperNodeIdBound(); ++v) {
        if (G.hasNode(v) &&
            dijkstra.getDistance(v) != std::numeric_limits<double>::infinity()) {
            totalDist += dijkstra.getDistance(v);
            reachableNodes++;
        }
    }

    std::cout << "Average distance: "
              << (reachableNodes > 0 ? totalDist / reachableNodes : 0) << std::endl;
    std::cout << "Reachable nodes: " << reachableNodes
              << " out of " << G.numberOfNodes() << std::endl;
}

void run_test(
        const std::string &label, NetworKit::Graph &G, const std::string &algorithm_name,
        Algorithm algorithm) {
    std::cout << "Graph: " << label << std::endl;
    std::cout << "Nodes: " << G.numberOfNodes()
              << ", Edges: " << G.numberOfEdges() << std::endl;

    runDijkstra(G, 0, algorithm_name, algorithm);
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <queue_type> <graph_file>" << std::endl;
        std::cerr << "Available queue types: ";
        for (const auto& [name, _] : ALGORITHM) {
            std::cerr << name << " ";
        }
        std::cerr << std::endl;
        return 1;
    }

    std::string algorithm_name(argv[1]);
    auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument("Unknown algorithm: " + algorithm_name);
    }
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph G) {
        run_test(label, G, algorithm_name, algorithm->second);
    });

    return 0;
}
