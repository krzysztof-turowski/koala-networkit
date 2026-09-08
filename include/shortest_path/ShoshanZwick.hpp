#pragma once

#include <utility>
#include <networkit/graph/Graph.hpp>
#include <Eigen/Dense>
#include "shortest_path/AllPairsShortestPaths.hpp"

namespace Koala {

/**
 * @ingroup shortest_path
 * Exact All-Pairs Shortest Paths for unweighted, undirected graphs
 * using the Shoshan-Zwick algorithm: repeated boolean matrix multiplication
 * with distance bit-recovery.
 *
 * Time complexity: O(M(n) * log n) where M(n) is the boolean
 * matrix multiplication time.
 */
class ShoshanZwickAPSP : public AllPairsShortestPaths<int> {
 public:
    using Matrix = Eigen::MatrixXi;

    explicit ShoshanZwickAPSP(const NetworKit::Graph &graph): AllPairsShortestPaths<int>(graph) { checkInput(); }

    explicit ShoshanZwickAPSP(NetworKit::Graph &&graph): AllPairsShortestPaths<int>(std::move(graph)) { checkInput(); }

    void run() override;

 private:
    // Throw if the stored graph is directed or weighted.
    void checkInput() const;

    static Matrix boolMul(const Matrix& A, const Matrix& B);
};

}  // namespace Koala
