#pragma once

#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup shortest_path
 * Approximate All-Pairs Shortest Paths for unweighted, undirected graphs
 * using the algorithm of Aingworth, Chekuri, Indyk, Motwani.
 *
 * Produces an additive +2 approximation: for every pair (u,v),
 *   d(u,v) <= d_hat(u,v) <= d(u,v) + 2
 *
 * Time complexity: O(n^{2.5} * sqrt(log n)) with threshold s = ceil(sqrt(n * log n)).
 */
class AingworthAlgorithm : public NetworKit::Algorithm {
 public:
    explicit AingworthAlgorithm(const NetworKit::Graph& G);

    void run() override;

    const std::vector<std::vector<int>>& getDistances() const;
    int getDistance(NetworKit::node u, NetworKit::node v) const;
    int getDiameter() const;

 private:
    const NetworKit::Graph* G;
    std::vector<std::vector<int>> distances;
    int diameter;

    static std::vector<int> bfs(
        NetworKit::node src, const std::vector<std::vector<NetworKit::node>>& adj,
        NetworKit::count n);

    static std::vector<int> bfsLow(
        NetworKit::node src, const std::vector<std::vector<NetworKit::node>>& adj,
        const std::vector<bool>& is_low, NetworKit::count n);

    static std::vector<NetworKit::node> dominatingSet(
        const std::vector<NetworKit::node>& high,
        const std::vector<std::vector<NetworKit::node>>& adj,
        const std::vector<bool>& is_high, NetworKit::count n);
};

}  // namespace Koala
