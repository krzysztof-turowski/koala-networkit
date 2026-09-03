#pragma once

#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup shortest_path
 * All-Pairs Shortest Paths for weighted, undirected graphs
 * using Spira's algorithm.
 *
 * Time complexity: O(n^2 * m * log n) worst case.
 */
class SpiraAlgorithm : public NetworKit::Algorithm {
 public:
    explicit SpiraAlgorithm(const NetworKit::Graph& G);

    void run() override;

    const std::vector<std::vector<NetworKit::edgeweight>>& getDistances() const;
    NetworKit::edgeweight getDistance(NetworKit::node u, NetworKit::node v) const;

 private:
    const NetworKit::Graph* G;
    std::vector<std::vector<NetworKit::edgeweight>> distances;
};

}  // namespace Koala
