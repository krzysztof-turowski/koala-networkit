#include <NTL/mat_ZZ_p.h>
#include <NTL/vector.h>

#include <iostream>

#include "matching/gaussian_matching/NaiveGaussElimination.hpp"
#include "matching/gaussian_matching/utils.hpp"

namespace Koala {
bool eliminate(MatZp &A, int r, int c) {
  if (A[r][c] == 0)
    return false;

  auto B = zeroMat(A.NumCols(), A.NumRows());
  for (int i = 0; i < A.NumCols(); i++) {
    for (int j = 0; j < A.NumRows(); j++) {
      B[i][j] = A[r][j] * A[i][c] / A[r][c];
    }
  }
  A -= B;
  return true;
}

std::vector<int> NaiveGaussElimination::pivotElimination(
    MatZp &A, std::function<bool(int, int)> isCellAllowed, bool bipartite = false) {
  int n = A.NumCols();

  std::vector<int> res(n);
  for (int c = 0; c < n; ++c) {
    for (int r = 0; r < n; ++r) {
      if (A[r][c] == 0 || !isCellAllowed(c, r))
        continue;

      eliminate(A, r, c);
      if (!bipartite) {
        eliminate(A, c, r);
      }

      res[c] = r;
      break;
    }
  }

  return res;
}

std::vector<int> NaiveGaussElimination::simpleElimination(MatZp &A, int k) {
  std::vector<int> res;

  for (int i = 0; i < k; ++i) {
    if (eliminate(A, i, i)) {
      res.push_back(i);
    }
  }

  return res;
}
}  // namespace Koala
