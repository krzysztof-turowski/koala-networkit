#pragma once

#include <utility>
#include <networkit/graph/Graph.hpp>
#include <Eigen/Dense>
#include "shortest_path/AllPairsShortestPaths.hpp"

namespace Koala {

/**
 * @ingroup shortest_path
 * All-Pairs Shortest Paths for unweighted, undirected graphs
 * using Seidel's APD (All Pairs Distances) algorithm.
 *
 * Time complexity: O(M(n) * log n) where M(n) is the matrix
 * multiplication time.
 */
class SeidelAPSP : public AllPairsShortestPaths<int> {
 public:
    using Matrix = Eigen::MatrixXi;

    explicit SeidelAPSP(const NetworKit::Graph &graph): AllPairsShortestPaths<int>(graph) { checkInput(); }

    explicit SeidelAPSP(NetworKit::Graph &&graph): AllPairsShortestPaths<int>(std::move(graph)) { checkInput(); }

    void run() override;

 private:
    // Throw if the stored graph is directed or weighted.
    void checkInput() const;

    Matrix APD(const Matrix& A) const;
};

}  // namespace Koala
