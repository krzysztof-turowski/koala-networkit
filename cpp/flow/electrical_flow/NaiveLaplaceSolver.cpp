#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "flow/electrical_flow/LaplaceSolver.hpp"

namespace Koala {
Eigen::VectorXd solveLaplace(
    const NetworKit::Graph &graph, const std::vector<std::vector<double>> &weights,
    const Eigen::VectorXd &b) {
  int N = graph.numberOfNodes();

  Eigen::MatrixXd L(N, N);
  for (int u = 0; u < N; ++u) {
    double wuu = 0;
    for (int v = 0; v < N; ++v) {
      if (graph.hasEdge(u, v)) {
        L(u, v) = -weights[u][v];
        wuu += weights[u][v];
      } else {
        L(u, v) = 0;
      }
    }
    L(u, u) = wuu;
  }

  auto Linv = L.completeOrthogonalDecomposition().pseudoInverse();
  return Linv * b;
}

}  // namespace Koala
