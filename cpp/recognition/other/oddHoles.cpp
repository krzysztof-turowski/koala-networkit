#include "oddHoles.h"
#include <algorithm>
#include <set>
#include "commons.h"
#include "testCommons.h"

using std::get;
using std::to_string;

vec<int> findHoleOfSize(const Graph &G, int size) {
  if (size <= 3) return vec<int>();

  vec<int> v;
  while (true) {
    nextPathInPlace(G, v, size, true, false, true);
    if (v.size() == size && isHole(G, v)) return v;

    if (v.empty()) break;
  }

  return vec<int>();
}

bool constainsHoleOfSize(const Graph &G, int size) { return !findHoleOfSize(G, size).empty(); }

vec<int> findOddHoleNaive(const Graph &G, bool gatherStats) {
  // for (int size = 5; size <= G.n; size += 2) {
  //   if (gatherStats) StatsFactory::startTestCasePart(to_string(size));
  //   auto v = findHoleOfSize(G, size);
  //   if (!v.empty()) return v;
  // }

  // return vec<int>();
  vec<int> v;
  nextPathInPlace(G, v, 0, true, false, true);
  return v;
}

bool containsOddHoleNaive(const Graph &G, bool gatherStats) {
  // auto v = findOddHoleNaive(G, gatherStats);
  // cout<<v<<endl;
  // return !v.empty();
  return !findOddHoleNaive(G, gatherStats).empty();
}

bool isT1(const Graph &G, const vec<int> &v) {
  if (v.size() != 5) return false;

  if (!isDistinctValues(v)) return false;

  for (int i : v)
    if (i < 0 || i >= G.n) return false;

  for (int i = 0; i < 5; i++) {
    for (int j = i + 1; j < 5; j++) {
      if (abs(i - j) == 1 || abs(i - j) == 4) {
        if (!G.areNeighbours(v[i], v[j])) return false;
      } else {
        if (G.areNeighbours(v[i], v[j])) return false;
      }
    }
  }

  return true;
}
vec<int> findT1(const Graph &G) { return findHoleOfSize(G, 5); }
bool containsT1(const Graph &G) { return !findT1(G).empty(); }

tuple<vec<int>, vec<int>, vec<int>> findT2(const Graph &G) {
  for (int v1 = 0; v1 < G.n; v1++) {
    for (int v2 : G[v1]) {
      for (int v3 : G[v2]) {
        if (v3 == v1) continue;
        for (int v4 : G[v3]) {
          if (v4 == v1 || v4 == v2) continue;
          if (!isAPath(G, vec<int>{v1, v2, v3, v4})) continue;
          auto Y = getCompleteVertices(G, {v1, v2, v4});
          auto antiCY = getComponentsOfInducedGraph(G.getComplement(), Y);

          for (auto X : antiCY) {
            auto P = findShortestPathWithPredicate(G, v1, v4, [&](int v) -> bool {
              if (v == v2 || v == v3) return false;
              if (G.areNeighbours(v, v2) || G.areNeighbours(v, v3)) return false;
              if (isComplete(G, X, v)) return false;

              return true;
            });

            if (!P.empty()) {
              return make_tuple(vec<int>{v1, v2, v3, v4}, P, X);
            }
          }
        }
      }
    }
  }

  return make_tuple(vec<int>(), vec<int>(), vec<int>());
}
bool containsT2(const Graph &G) { return !get<0>(findT2(G)).empty(); }

tuple<vec<int>, vec<int>, vec<int>> findT3(const Graph &G) {
  for (int v1 = 0; v1 < G.n; v1++) {
    for (int v2 : G[v1]) {
      for (int v5 = 0; v5 < G.n; v5++) {
        if (v5 == v1 || v5 == v2 || G.areNeighbours(v5, v1) || G.areNeighbours(v5, v2)) continue;

        auto Y = getCompleteVertices(G, {v1, v2, v5});
        auto antiCY = getComponentsOfInducedGraph(G.getComplement(), Y);
        for (auto X : antiCY) {
          set<int> Fprim;
          vec<int> visited(G.n);
          dfsWith(G, visited, v5, [&](int v) -> void { Fprim.insert(v); },
                  [&](int v) -> bool {
                    if (G.areNeighbours(v1, v) || G.areNeighbours(v2, v)) return false;
                    if (isComplete(G, X, v)) return false;

                    return true;
                  });

          set<int> F(Fprim.begin(), Fprim.end());
          for (int fp : Fprim) {
            for (int v : G[fp]) {
              if (F.count(v) == 0 && isComplete(G, X, v) && !G.areNeighbours(v, v1) &&
                  !G.areNeighbours(v, v2) && !G.areNeighbours(v, v5))
                F.insert(v);
            }
          }

          for (int v4 : G[v1]) {
            if (G.areNeighbours(v4, v2) || G.areNeighbours(v4, v5)) continue;
            int v6 = -1;
            for (int f : F) {
              if (G.areNeighbours(v4, f)) {
                v6 = f;
                break;
              }
            }
            if (v6 == -1) continue;

            bool v4HasNonNeighbourInX = false;
            for (int x : X) {
              if (!G.areNeighbours(v4, x)) {
                v4HasNonNeighbourInX = true;
                break;
              }
            }
            if (!v4HasNonNeighbourInX) continue;

            for (int v3 : G[v2]) {
              if (!G.areNeighbours(v3, v4) || !G.areNeighbours(v3, v5) || G.areNeighbours(v3, v1)) continue;

              bool v3HasNonNeighbourInX = false;
              for (int x : X) {
                if (!G.areNeighbours(v3, x)) {
                  v3HasNonNeighbourInX = true;
                  break;
                }
              }
              if (!v3HasNonNeighbourInX) continue;

              auto P =
                  findShortestPathWithPredicate(G, v6, v5, [&](int v) -> bool { return Fprim.count(v) > 0; });

              if (P.empty()) {  // Should not happen
                throw std::logic_error("Algorithm Error: Could not find path P in T3.");
              }

              return make_tuple(vec<int>{v1, v2, v3, v4, v5, v6}, P, X);
            }
          }
        }
      }
    }
  }

  return make_tuple(vec<int>(), vec<int>(), vec<int>());
}
bool containsT3(const Graph &G) { return !get<0>(findT3(G)).empty(); }

bool isT3(const Graph &G, const vec<int> &v_const, const vec<int> &P, const vec<int> &X_const) {
  if (v_const.size() != 6 || P.empty() || X_const.empty()) return false;

  vec<int> v(v_const.begin(), v_const.end());
  v.insert(v.begin(), -1);  // To have v1,... v6

  if (!isDistinctValues(v)) return false;

  auto edges = vec<vec<int>>{{1, 2}, {3, 4}, {1, 4}, {2, 3}, {3, 5}, {4, 6}};
  auto nonEdges = vec<vec<int>>{{1, 3}, {2, 4}, {1, 5}, {2, 5}, {1, 6}, {2, 6}, {4, 5}};

  for (auto e : edges) {
    if (!G.areNeighbours(v[e[0]], v[e[1]])) return false;
  }

  for (auto e : nonEdges) {
    if (G.areNeighbours(v[e[0]], v[e[1]])) return false;
  }

  auto antiC = getComponentsOfInducedGraph(G.getComplement(), getCompleteVertices(G, {v[1], v[2], v[5]}));

  vec<int> X(X_const.begin(), X_const.end());
  sort(X.begin(), X.end());

  bool isXAnticomponent = false;
  for (auto ac : antiC) {
    sort(ac.begin(), ac.end());
    if (X == ac) {
      isXAnticomponent = true;
      break;
    }
  }
  if (!isXAnticomponent) return false;

  if (isComplete(G, X, v[3]) || isComplete(G, X, v[4])) return false;

  if ((P[0] != v[5] || P.back() != v[6]) && (P[0] != v[6] || P.back() != v[5])) return false;

  if (!isAPath(G, P)) return false;

  for (int i = 1; i < P.size() - 1; i++) {
    int u = P[i];
    for (int i = 1; i <= 4; i++) {
      if (u == v[i]) return false;
    }
    for (int x : X) {
      if (u == x) return false;
    }

    if (isComplete(G, X, u)) return false;

    if (G.areNeighbours(v[1], u) || G.areNeighbours(v[2], u)) return false;
  }

  if (G.areNeighbours(v[5], v[6]) && isComplete(G, X, v[6])) return false;

  return true;
}
