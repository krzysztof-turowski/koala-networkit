#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <Eigen/Dense>

namespace Koala {

/**
 * @ingroup shortest_path
 * All-Pairs Shortest Paths for unweighted, undirected graphs
 * using Seidel's APD (All Pairs Distances) algorithm.
 *
 * Time complexity: O(M(n) * log n) where M(n) is the matrix
 * multiplication time.
 */
class APDAlgorithm : public NetworKit::Algorithm {
 public:
    using Matrix = Eigen::MatrixXi;

    explicit APDAlgorithm(const NetworKit::Graph& G);

    void run() override;

    const Matrix& getDistances() const;
    int getDistance(NetworKit::node u, NetworKit::node v) const;

 private:
    const NetworKit::Graph* G;
    Matrix distances;

    Matrix APD(const Matrix& A) const;
};

}  // namespace Koala
