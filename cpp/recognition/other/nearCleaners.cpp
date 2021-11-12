
#include <boost/dynamic_bitset.hpp>
#include <set>

#include "commons.h"
#include "nearCleaners.h"
#include "testCommons.h"

bool containsOddHoleWithNearCleanerX(const Graph &G, const set<int> &sX, const vec<vec<int>> &triplePaths,
                                     bool gatherStats) {
  if (gatherStats) StatsFactory::startTestCasePart("Test NC Shortest Paths");
  vec<vec<int>> penultimate;
  auto R = allShortestPathsWithPredicate(G, [&](int v) -> bool { return sX.count(v) == 0; }, penultimate);

  if (gatherStats) StatsFactory::startTestCasePart("Test NC Rest");
  for (int y1 = 0; y1 < G.n; y1++) {
    if (sX.count(y1) > 0) continue;

    for (auto x : triplePaths) {
      bool containsY1 = false;
      for (int i = 0; i < 3; i++) {
        if (x[i] == y1) {
          containsY1 = true;
          break;
        }
      }
      if (containsY1) continue;

      int x1 = x[0], x3 = x[1], x2 = x[2];
      if (R[x1][y1] == 0 || R[x2][y1] == 0) continue;

      int y2 = penultimate[x2][y1];

      int n = R[x2][y1];
      if (R[x1][y1] + 1 != n || R[x1][y2] != n) continue;

      if ((R[x3][y1] < n) || (R[x3][y2] < n)) continue;

      // TODO(Adrian) remove, use some code coverage tool instead
      // cout << "Interesting odd hole: " << endl;
      // cout << "sX: " << sX << endl;
      // cout << "x1: " << x1 << endl;
      // cout << "x2: " << x2 << endl;
      // cout << "x3: " << x3 << endl;
      // cout << "y1: " << y1 << endl;
      // cout << "y2: " << y2 << endl;
      // cout << "R(x1, y1): " << R[x1][y1] << endl;
      // cout << "R(x2, y1): " << R[x2][y1] << endl;
      // cout << "R(x3, y1): " << R[x3][y1] << endl;
      // cout << "R(x3, y2): " << R[x3][y2] << endl;
      // cout << G << endl;

      return true;
    }
  }

  return false;
}

bool isRelevantTriple(const Graph &G, int a, int b, int c) {
  // for (int i = 0; i < 3; i++) {
  //   if (v[i] < 0 || v[i] >= G.n) return false;
  // }

  // int a = v[0], b = v[1], c = v[2];

  if (a == b || G.areNeighbours(a, b)) return false;

  if (G.areNeighbours(a, c) && G.areNeighbours(b, c)) return false;

  return true;
}

boost::dynamic_bitset<ul> getXforRelevantTriple(const Graph &G, const Graph &GC, int a, int b, int c) {
  auto antiCompsNab = getComponentsOfInducedGraph(GC, getCompleteVertices(G, {a, b}));
  int r = 0;
  for (auto comp : antiCompsNab) {
    if (comp.size() <= r) continue;

    bool containsNonNofC = false;
    for (int v : comp) {
      if (!G.areNeighbours(c, v)) containsNonNofC = true;
    }

    if (containsNonNofC) r = comp.size();
  }

  vec<int> Y;
  for (auto comp : antiCompsNab) {
    if (comp.size() > r) {
      Y.insert(Y.end(), comp.begin(), comp.end());
    }
  }

  vec<int> W;
  for (auto comp : antiCompsNab) {
    for (int v : comp) {
      if (!G.areNeighbours(v, c)) {
        W = comp;
        break;
      }
    }
  }
  W.push_back(c);

  W.insert(W.end(), Y.begin(), Y.end());
  auto Z = getCompleteVertices(G, W);

  boost::dynamic_bitset<ul> ret = getBitset(G.n, Y);
  ret |= getBitset(G.n, Z);
  return ret;
}

set<boost::dynamic_bitset<ul>> getPossibleNearCleaners(const Graph &G, const Graph &GC, bool gatherStats) {
  if (gatherStats) StatsFactory::startTestCasePart("NC 1");

  vec<boost::dynamic_bitset<ul>> Ns;
  for (int u = 0; u < G.n; u++) {
    for (int v : G[u]) {
      Ns.push_back(getBitset(G.n, getCompleteVertices(G, {u, v})));
    }
  }

  if (gatherStats) StatsFactory::startTestCasePart("NC 2");

  vec<boost::dynamic_bitset<ul>> Xs;
  for (int a = 0; a < G.n; a++) {
    for (int b = 0; b < G.n; b++) {
      if (a == b || G.areNeighbours(a, b)) continue;
      for (int c = 0; c < G.n; c++) {
        if (isRelevantTriple(G, a, b, c)) {
          Xs.push_back(getXforRelevantTriple(G, GC, a, b, c));
        }
      }
    }
  }

  if (gatherStats) StatsFactory::startTestCasePart("NC 3");

  set<boost::dynamic_bitset<ul>> res;
  for (auto N : Ns) {
    for (auto X : Xs) {
      boost::dynamic_bitset<ul> tmp = X;
      tmp |= N;
      res.insert(tmp);
    }
  }

  return res;
}
