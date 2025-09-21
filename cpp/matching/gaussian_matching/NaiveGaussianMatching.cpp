#include <matching/gaussian_matching/NaiveGaussElimination.hpp>
#include <matching/gaussian_matching/NaiveGaussianMatching.hpp>
#include <matching/gaussian_matching/utils.hpp>

#include <NTL/ZZ_p.h>
#include <networkit/graph/Graph.hpp>

using namespace std;
using namespace NetworKit;

namespace Koala {
static MatZp generateMatrix(const Graph &G);

NaiveGaussianMatching::NaiveGaussianMatching(const Graph &G1) : G(G1) {}

Matching NaiveGaussianMatching::getMatching() { return M; }

void NaiveGaussianMatching::run() {
  M.clear();
  MatZp AG = generateMatrix(G);
  auto M1 = NaiveGaussElimination::pivotElimination(
      AG, [&AG](int r, int c) { return AG[r][c] != 0; });
  for (int c = 0; c < M1.size(); c++) {
    M.insert({c, M1[c]});
  }
}

static MatZp generateMatrix(const Graph &G) {
  int n = G.numberOfNodes();

  auto AG = zeroMat(n, n);
  for (auto [u, v] : G.edgeRange()) {
    auto Xuv = generateRandom();
    AG[u][v] = Xuv;
    AG[v][u] = -Xuv;
  }

  MatZp Ainv;
  inv(Ainv, AG);
  return Ainv;
}
} // namespace Koala