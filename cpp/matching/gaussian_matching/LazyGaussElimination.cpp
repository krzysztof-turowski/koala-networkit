#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vector.h>
#include <iostream>
#include <matching/gaussian_matching/LazyGaussElimination.hpp>
#include <matching/gaussian_matching/utils.hpp>
#include <vector>

using namespace std;
using namespace NTL;

namespace Koala {
static pair<int, int> get2div(int n) {
  // return m such that largest m: 2^m|n
  int m = 0, m2 = 1;
  while (n % (m2 * 2) == 0)
    m = m + 1, m2 *= 2;
  return {m, m2};
}

static MatZp update(int r1, int r2, int c1, int c2,
                    vector<tuple<VecZp, VecZp, ZZ_p>> acc) {
  int k = acc.size();
  MatZp U, V, R;
  U.SetDims(r2 - r1 + 1, k);
  V.SetDims(k, c2 - c1 + 1);
  for (int i = 0; i < k; ++i) {
    auto [u, v, a] = acc[i];
    setCol(U, i, divVec(segment(u, r1, r2 + 1), a));
    setRow(V, i, segment(v, c1, c2 + 1));
  }
  return U * V;
}

vector<int>
LazyGaussElimination::pivotElimination(MatZp &A,
                                       function<bool(int, int)> isCellAllowed) {
  int n = A.NumCols();

  vector<tuple<VecZp, ZZ_p>> lazy(n);
  vector<int> order(n);
  for (int c = 0; c < n; ++c) {
    int r;
    for (r = 0; r < n; ++r) {
      if (A[r][c] != 0 && isCellAllowed(c, r))
        break;
    }
    assert(r != n);
    order[c] = r;

    auto [j, j2] = get2div(c + 1);

    auto a = A[r][c];
    lazy[c] = {getCol(A, c), a};

    vector<tuple<VecZp, VecZp, ZZ_p>> superLazy(j2);
    for (int i = 0; i < j2; ++i) {
      auto [l, l2] = get2div(i + 1);

      auto from = superLazy.begin() + max(0, i - l2 + 1);
      auto to = superLazy.begin() + i;
      auto acc = vector<tuple<VecZp, VecZp, ZZ_p>>(from, to);

      int r1 = 0, r2 = n - 1;
      int c1 = c + 1, c2 = min(c + j2, n - 1);

      if (r1 <= r2 && c1 <= c2) {
        auto B = update(0, n - 1, c1, c2, acc);
        subBlock(A, r1, c1, B);
      }

      auto v1 = getRow(A, order[c - j2 + 1 + i]);
      auto [u1, a1] = lazy[c - j2 + 1 + i];
      superLazy[i] = {u1, v1, a1};
    }

    VecZp zeroVec;
    zeroVec.SetLength(A.NumRows());
    setCol(A, c, zeroVec);
    A[r][c] = a;
  }
  return order;
}

vector<int> LazyGaussElimination::simpleElimination(MatZp &A, int k) {
  int n = A.NumCols();
  vector<int> res;

  vector<tuple<VecZp, VecZp, ZZ_p>> lazy(k);
  for (int i = 0; i < k; ++i) {
    if (A[i][i] == 0) {
      lazy[i] = {zeroVec(A.NumCols()), zeroVec(A.NumRows()), ZZ_p(1)};
    } else {
      lazy[i] = {getCol(A, i), getRow(A, i), A[i][i]};
      res.push_back(i);
    }
    auto a = A[i][i];

    auto [j, j2] = get2div(i + 1);

    auto from = lazy.begin() + max(0, i - j2 + 1);
    auto to = lazy.begin() + i + 1;
    auto acc = vector<tuple<VecZp, VecZp, ZZ_p>>(from, to);

    int r1 = i + 1, r2 = min(i + j2, n - 1);
    int c1 = i + 1, c2 = n - 1;
    if (r1 <= r2 && c1 <= c2) {
      auto B = update(r1, r2, c1, c2, acc);
      subBlock(A, r1, c1, B);
    }

    r1 = i + j2 + 1, r2 = n - 1;
    c1 = i + 1, c2 = min(i + j2, n - 1);
    if (r1 <= r2 && c1 <= c2) {
      auto B = update(r1, r2, c1, c2, acc);
      subBlock(A, r1, c1, B);
    }

    setCol(A, i, zeroVec(A.NumRows()));
    setRow(A, i, zeroVec(A.NumCols()));
    A[i][i] = a;
  }

  return res;
}
}  // namespace Koala
