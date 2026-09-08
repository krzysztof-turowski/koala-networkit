#pragma once

#include <utility>
#include <vector>
#include <networkit/graph/Graph.hpp>
#include "shortest_path/AllPairsShortestPaths.hpp"

namespace Koala {

/**
 * @ingroup shortest_path
 * All-Pairs Shortest Paths for weighted, undirected graphs
 * using Spira's algorithm (run once per source).
 *
 * Time complexity: O(n * m * log n) worst case, O(n^2 * log n) average case.
 */
class SpiraAPSP : public AllPairsShortestPaths<NetworKit::edgeweight> {
 public:
    explicit SpiraAPSP(const NetworKit::Graph &graph): AllPairsShortestPaths<NetworKit::edgeweight>(graph) { checkInput(); }

    explicit SpiraAPSP(NetworKit::Graph &&graph): AllPairsShortestPaths<NetworKit::edgeweight>(std::move(graph)) { checkInput(); }

    void run() override;

 private:
    // Throw if the stored graph is not weighted.
    void checkInput() const;
};

}  // namespace Koala
