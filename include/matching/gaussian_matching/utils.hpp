#pragma once

#include "networkit/graph/GraphTools.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <networkit/graph/Graph.hpp>
namespace Koala {
constexpr double EPS = 1e-8;
constexpr int ZP_MOD = 7727;

typedef NTL::ZZ_p Zp;
typedef NTL::Vec<Zp> VecZp;
typedef NTL::Mat<Zp> MatZp;

inline bool eq(double a, double b) { return fabs(a - b) <= EPS; }

inline void initZp(int p) { Zp::init(NTL::conv<NTL::ZZ>(p)); }

inline Zp generateRandom() {
  unsigned int seed = (unsigned int)time(NULL);
  return Zp(rand_r(&seed));
}

inline VecZp getCol(const MatZp &A, int c) {
  VecZp col;
  col.SetLength(A.NumRows());
  for (int i = 0; i < col.length(); ++i) {
    col[i] = A[i][c];
  }
  return col;
}

inline void setCol(MatZp &A, int c, const VecZp &vec) {
  for (int i = 0; i < vec.length(); ++i) {
    A[i][c] = vec[i];
  }
}

inline VecZp getRow(const MatZp &A, int r) {
  VecZp row;
  row.SetLength(A.NumCols());
  for (int i = 0; i < row.length(); ++i) {
    row[i] = A[r][i];
  }
  return row;
}

inline void setRow(MatZp &A, int r, const VecZp &vec) {
  for (int i = 0; i < vec.length(); ++i) {
    A[r][i] = vec[i];
  }
}

inline VecZp divVec(const VecZp &vec, Zp a) {
  VecZp res;
  res.SetLength(vec.length());
  for (int i = 0; i < vec.length(); ++i) {
    res[i] = vec[i] / a;
  }
  return res;
}

inline VecZp segment(const VecZp &vec, int b, int e) {
  VecZp seg;
  seg.SetLength(e - b);
  for (int i = b; i < e; ++i) {
    seg[i - b] = vec[i];
  }
  return seg;
}

inline VecZp zeroVec(int n) {
  VecZp vec;
  vec.SetLength(n);
  return vec;
}

inline MatZp zeroMat(int r, int c) {
  MatZp A;
  A.SetDims(r, c);
  return A;
}

inline void subBlock(MatZp &A, int r1, int c1, const MatZp &B) {
  for (int i = 0; i < B.NumRows(); ++i) {
    for (int j = 0; j < B.NumCols(); ++j) {
      A[r1 + i][c1 + j] -= B[i][j];
    }
  }
}

inline std::pair<NetworKit::Graph, std::vector<int>>
reindexGraph(const NetworKit::Graph &G) {
  auto indexes = NetworKit::GraphTools::getContinuousNodeIds(G);
  auto G1 = NetworKit::GraphTools::getCompactedGraph(G, indexes);

  std::vector<int> labels(G1.numberOfNodes());
  for (auto [k, v] : indexes) {
    labels[v] = k;
  }

  return {G1, labels};
}

inline std::tuple<NetworKit::Graph, std::vector<int>, MatZp>
reindexGraph(const NetworKit::Graph &G, const MatZp &AG) {
  auto [G1, labels] = reindexGraph(G);

  int n = G.numberOfNodes();
  auto AG1 = zeroMat(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      AG1[i][j] = AG[labels[i]][labels[j]];
    }
  }

  return {G1, labels, AG1};
}
}  // namespace Koala
