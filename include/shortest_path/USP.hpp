#pragma once

#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <Eigen/Dense>

namespace Koala {

/**
 * @ingroup shortest_path
 * All-Pairs Shortest Paths for unweighted, undirected graphs
 * using the USP (Unweighted Shortest Paths) algorithm based on
 * repeated boolean matrix multiplication.
 *
 * Time complexity: O(M(n) * log n) where M(n) is the boolean
 * matrix multiplication time.
 */
class USPAlgorithm : public NetworKit::Algorithm {
 public:
    using Matrix = Eigen::MatrixXi;

    explicit USPAlgorithm(const NetworKit::Graph& G);

    void run() override;

    const Matrix& getDistances() const;
    int getDistance(NetworKit::node u, NetworKit::node v) const;

 private:
    const NetworKit::Graph* G;
    Matrix distances;

    static Matrix boolMul(const Matrix& A, const Matrix& B);
};

}  // namespace Koala
