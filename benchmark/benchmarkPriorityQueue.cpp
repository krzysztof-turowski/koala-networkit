#include <cassert>
#include <iostream>
#include <map>
#include <chrono>
#include <limits>
#include <vector>
#include <string>
#include <queue>
#include <iomanip>
#include <type_traits>

#include <networkit/graph/Graph.hpp>
#include <io/DimacsGraphReader.hpp>

#include "structures/PriorityQueue.hpp"

#include "structures/priority_queue/RankPairingHeap.hpp"
#include "structures/priority_queue/SkewHeap.hpp"
#include "structures/priority_queue/WeakHeap.hpp"

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
            NetworKit::node u = pq.pop();

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

std::map<std::string, int> QUEUE_TYPE = {
    { "rank_pairing_heap", 0 },
    { "skew_heap", 1 },
    { "weak_heap", 2 }
};

void runDijkstra(NetworKit::Graph& G, NetworKit::node source, const std::string& queueType) {
    auto start = std::chrono::high_resolution_clock::now();

    Dijkstra dijkstra(G, source);

    auto it = QUEUE_TYPE.find(queueType);
    if (it == QUEUE_TYPE.end()) {
        std::cerr << "Unknown queue type: " << queueType << std::endl;
        return;
    }

    switch (it->second) {
        case 0:
            dijkstra.run<Koala::RankPairingHeap<NetworKit::node>>();
            break;
        case 1:
            dijkstra.run<Koala::SkewHeap<NetworKit::node>>();
            break;
        case 2:
            dijkstra.run<Koala::WeakHeap<NetworKit::node>>();
            break;
        default:
            std::cerr << "Unhandled queue type: " << queueType << std::endl;
            return;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Algorithm: Dijkstra with " << queueType << std::endl;
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

void run_dimacs_test(const std::string &path, const std::string &queueType) {
    auto G = Koala::DimacsGraphReader().read(path);

    std::cout << "Graph: " << path << std::endl;
    std::cout << "Nodes: " << G.numberOfNodes()
              << ", Edges: " << G.numberOfEdges() << std::endl;

    runDijkstra(G, 0, queueType);
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <queue_type> <graph_file>" << std::endl;
        std::cerr << "Available queue types: ";
        for (const auto& [name, _] : QUEUE_TYPE) {
            std::cerr << name << " ";
        }
        std::cerr << std::endl;
        return 1;
    }

    std::string queueType(argv[1]);
    std::string path(argv[2]);

    run_dimacs_test(path, queueType);

    return 0;
}
