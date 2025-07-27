#include <Eigen/Core>
#include <Eigen/Dense>
#include <flow/electric_flow/LaplaceSolver.hpp>

using namespace Eigen;

namespace Koala {
VectorXd solveLaplace(const Graph &graph, const vector<vector<double>> &weights,
                      const VectorXd &b) {
  int N = graph.numberOfNodes();

  MatrixXd L(N, N);
  for (int u = 0; u < N; ++u) {
    double wuu=0;
    for (int v = 0; v < N; ++v) {
      if (graph.hasEdge(u, v)) {
        L(u, v) = -weights[u][v];
        wuu+=weights[u][v];
      } else {
        L(u, v) = 0;
      }
    }
    L(u,u) = wuu;
  }

  auto Linv = L.completeOrthogonalDecomposition().pseudoInverse();
  VectorXd x = Linv * b;
  // cout << L << endl;
  // cout << b << endl;
  // cout << x << endl;
  return x;
}

} // namespace Koala
