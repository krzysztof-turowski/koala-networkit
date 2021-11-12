#include "pyramids.h"
#include <functional>
#include <queue>
#include <unordered_set>
#include "commons.h"

using std::get;
using std::unordered_set;

bool checkPrerequisites(const Graph &G, const vec<int> &b, const int a, const vec<int> &s) {
  // We assume b is a Triangle and s is an EmptyStarTriangle

  int aAdjB = 0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) continue;
      if (b[i] == s[j]) {
        return false;
      }
      if ((G.areNeighbours(b[i], s[j]) && b[j] != s[j]) || G.areNeighbours(s[i], s[j])) {
        return false;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    if (!G.areNeighbours(a, s[i])) return false;

    if (G.areNeighbours(a, b[i])) {
      aAdjB++;
      if (aAdjB > 1) return false;
      if (b[i] != s[i]) return false;
    }
  }

  return true;
}

bool vectorsCutEmpty(vec<int>::iterator aBegin, vec<int>::iterator aEnd, vec<int>::iterator bBegin,
                     vec<int>::iterator bEnd) {
  unordered_set<int> a;
  a.reserve(aEnd - aBegin);
  for (auto it = aBegin; it != aEnd; it++) {
    a.insert(*it);
  }

  for (auto it = bBegin; it != bEnd; it++) {
    if (a.count(*it) > 0) return false;
  }

  return true;
}
bool noEdgesBetweenVectors(const Graph &G, vec<int>::iterator aBegin, vec<int>::iterator aEnd,
                           vec<int>::iterator bBegin, vec<int>::iterator bEnd) {
  unordered_set<int> aTouches;
  aTouches.reserve(aEnd - aBegin);
  for (auto it = aBegin; it != aEnd; it++) {
    for (int i : G[*it]) {
      aTouches.insert(i);
    }
  }

  for (auto it = bBegin; it != bEnd; it++) {
    if (aTouches.count(*it) > 0) return false;
  }

  return true;
}

bool containsPyramid(const Graph &G) {
  auto t = findPyramid(G);
  return get<0>(t).size() > 0;
}

vec<vec<vec<int>>> calculatePPaths(const Graph &G, const vec<int> &b, const vec<int> &s, const vec<int> &M) {
  vec<vec<int>> S[3];  // S[a][b][c] = S_a(b)[c], c-th vertex of the Sa(b) path
  vec<vec<int>> T[3];
  for (int i = 0; i < 3; i++) {
    S[i].resize(G.n);
    T[i].resize(G.n);

    for (int m = 0; m < G.n; m++) {
      // if (!M[m])
      //   continue;

      auto noNeighbours = [&](int v) {
        if (v == s[i] || v == m) return true;
        if (!M[v]) return false;
        for (int j = 0; j < 3; j++) {
          if (i == j) continue;
          // if (b[i] == v || s[i] == v)
          // continue;
          if (G.areNeighbours(b[j], v) || G.areNeighbours(s[j], v)) return false;
        }

        return true;
      };

      S[i][m] = findShortestPathWithPredicate(G, s[i], m, noNeighbours);
      T[i][m] = findShortestPathWithPredicate(G, m, b[i], noNeighbours);
      // TODO(Adrian) test it somehow
    }
  }

  // cout << "s: " << s << endl;
  // cout << "b: " << b << endl;

  vec<vec<vec<int>>> P(3);
  for (int i = 0; i < 3; i++) {
    P[i].resize(G.n);
    if (s[i] == b[i]) {
      P[i][b[i]] = vec<int>{b[i]};
    } else {
      for (int m = 0; m < G.n; m++) {
        // if (!M[m]) {
        // cout << endl;
        // continue;
        // }

        bool pathExists = true;

        // Check if S and T paths exist or are between same vertices
        if ((S[i][m].size() == 0 && s[i] != m) || (T[i][m].size() == 0 && m != b[i])) {
          pathExists = false;
        }

        for (int j = 0; pathExists && j < 3; j++) {
          if (i == j) continue;

          if (M[m] && (G.areNeighbours(m, s[j]) || G.areNeighbours(m, b[j]))) {
            pathExists = false;
          }
        }

        if (pathExists &&
            !vectorsCutEmpty(S[i][m].begin(), S[i][m].end() - 1, T[i][m].begin() + 1, T[i][m].end())) {
          pathExists = false;
        }

        if (pathExists &&
            !noEdgesBetweenVectors(G, S[i][m].begin(), S[i][m].end() - 1, T[i][m].begin() + 1,
                                   T[i][m].end())) {
          pathExists = false;
        }

        if (pathExists) {
          P[i][m].insert(P[i][m].end(), S[i][m].begin(), S[i][m].end());
          P[i][m].insert(P[i][m].end(), T[i][m].begin() + 1, T[i][m].end());
        }
        // if (P[i][m].size() > 0 || S[i][m].size() > 0 || T[i][m].size() > 0) {
        //   cout << "S_" << i << "(" << s[i] << "," << m << ")= " << S[i][m] << "\t";
        //   cout << "T_" << i << "(" << m << "," << b[i] << ")= " << T[i][m] << "\t";
        //   cout << "  P[" << i << "][" << m << "](" << s[i] << "," << b[i] << ")= " << P[i][m] << endl;
        // }
      }
    }
  }  // P calculated
  return P;
}

tuple<vec<int>, int, vec<vec<int>>> findPyramid(const Graph &G) {
  auto triangles = getTriangles(G);
  auto emptyStars = getEmptyStarTriangles(G);

  for (auto triangle : triangles) {
    for (auto eStar : emptyStars) {
      const vec<int> &b = triangle;
      const int a = eStar.st;
      const vec<int> &s = eStar.nd;

      if (!checkPrerequisites(G, b, a, s)) continue;

      vec<int> M(G.n, 1);
      for (int i = 0; i < 3; i++) M[s[i]] = M[b[i]] = 0;

      vec<vec<vec<int>>> P = calculatePPaths(G, b, s, M);

      vec<vec<int>> pairs{{0, 1}, {0, 2}, {1, 2}};
      vec<vec<int>> goodPairs[3];
      for (int k = 0; k < 3; k++) {
        goodPairs[k].resize(G.n);
        for (int i = 0; i < G.n; i++) goodPairs[k][i].resize(G.n);
      }

      vec<int> color(G.n);

      for (int iPair = 0; iPair < 3; iPair++) {
        auto &pair = pairs[iPair];
        for (int m1 = 0; m1 < G.n; m1++) {
          if (P[pair[0]][m1].size() == 0) continue;

          for (int i = 0; i < G.n; i++) color[i] = 0;  // 0 - white

          for (int i : P[pair[0]][m1]) {
            if (!M[i]) continue;
            color[i] = 1;
            for (int j : G[i]) color[j] = 1;
          }

          for (int m2 = 0; m2 < G.n; m2++) {
            if (P[pair[1]][m2].size() == 0) continue;

            bool isOk = true;
            for (int i : P[pair[1]][m2])
              if (color[i]) {
                isOk = false;
                break;
              }

            if (isOk) {
              goodPairs[iPair][m1][m2] = true;
            }
          }
        }
      }  // (i, j) good pairs completed

      vec<vec<int>> triples = generateTuples(3, G.n);
      for (auto &triple : triples) {
        bool isOk = true;
        for (int iPair = 0; iPair < 3; iPair++) {
          auto &pair = pairs[iPair];
          if (!goodPairs[iPair][triple[pair[0]]][triple[pair[1]]]) {
            isOk = false;
            break;
          }
        }

        if (isOk) {
          return make_tuple(b, a, vec<vec<int>>{P[0][triple[0]], P[1][triple[1]], P[2][triple[2]]});
        }
      }
    }
  }

  return make_tuple(vec<int>(), -1, vec<vec<int>>());
}

bool isPyramid(const Graph &G, vec<int> b, int a, vec<vec<int>> P) {
  if (b.size() != 3 || a < 0 || a >= G.n || P.size() != 3) return false;

  for (int k = 0; k < 3; k++) {
    if (P[k].size() == 0) return false;
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) continue;
      if (!G.areNeighbours(b[i], b[j])) return false;
      if (b[i] == b[j]) return false;
      if (P[i][0] == P[j][0]) return false;
      if (G.areNeighbours(P[i][0], P[j][0])) return false;
    }
  }

  for (int i = 0; i < 3; i++) {
    if (!G.areNeighbours(a, P[i][0])) return false;

    if (P[i].back() != b[i]) return false;
  }

  for (int k = 0; k < 3; k++) {
    for (int i = 1; i < P[k].size(); i++) {
      if (!G.areNeighbours(P[k][i - 1], P[k][i])) return false;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = i + 1; j < 3; j++) {
      int numEdges = 0;
      for (int a : P[i]) {
        for (int b : P[j]) {
          if (G.areNeighbours(a, b)) numEdges++;
        }
      }

      if (numEdges != 1) return false;
    }
  }

  int adjA = 0;
  for (int i = 0; i < 3; i++) {
    if (G.areNeighbours(a, b[i])) adjA++;
  }
  if (adjA > 1) return false;

  return true;
}
