#pragma once

#include <vector>
#include <networkit/graph/Graph.hpp>
#include "shortest_path/AllPairsShortestPaths.hpp"

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
class AingworthAPSP : public AllPairsShortestPaths<int> {
 public:
    explicit AingworthAPSP(const NetworKit::Graph &graph): AllPairsShortestPaths<int>(graph) { checkInput(); }

    explicit AingworthAPSP(NetworKit::Graph &&graph): AllPairsShortestPaths<int>(std::move(graph)) { checkInput(); }

    /**
     * Use a custom degree threshold @p s (splitting low- and high-degree nodes)
     * instead of the default ceil(sqrt(n * log n)). Any value is accepted; a
     * non-positive @p s falls back to the default threshold.
     */
    AingworthAPSP(const NetworKit::Graph &graph, int s): AllPairsShortestPaths<int>(graph), threshold(s) { checkInput(); }

    AingworthAPSP(NetworKit::Graph &&graph, int s): AllPairsShortestPaths<int>(std::move(graph)), threshold(s) { checkInput(); }

    void run() override;

 private:
    // Degree threshold s; when <= 0 the default ceil(sqrt(n * log n)) is used.
    int threshold = 0;

    // Throw if the stored graph is directed or weighted.
    void checkInput() const;

    static std::vector<int> bfs(NetworKit::node src, const NetworKit::Graph& graph, NetworKit::count n);

    static std::vector<NetworKit::node> dominatingSet(
        const std::vector<std::vector<NetworKit::node>>& adj,
        const std::vector<bool>& is_high, NetworKit::count n);
};

}  // namespace Koala
