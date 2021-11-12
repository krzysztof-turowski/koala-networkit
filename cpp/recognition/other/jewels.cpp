#include "jewels.h"
#include <set>
#include "commons.h"

bool isJewel(const Graph &G, const vec<int> &v) {
  if (v.size() != 5) return false;

  set<int> S;
  for (int i : v) S.insert(i);
  if (S.size() != 5) return false;

  for (int i = 0; i < 5; i++) {
    if (!G.areNeighbours(v[i], v[(i + 1) % 5])) return false;
  }

  if (G.areNeighbours(v[0], v[2]) || G.areNeighbours(v[1], v[3]) || G.areNeighbours(v[0], v[3])) return false;

  auto noNeighbours = [&](int p) {
    if (p == v[1] || p == v[2] || p == v[4]) return false;

    for (int i : G[p]) {
      if (i == v[1] || i == v[2] || i == v[4]) return false;
    }
    return true;
  };

  vec<int> P = findShortestPathWithPredicate(G, v[0], v[3], noNeighbours);

  if (P.empty()) return false;

  // std::cout << "Jewel: ";
  // for (auto v : P)
      // std::cout << v << " ";
  // std::cout << std::endl;
  return true;
}

// returns [v1, ..., v5] or empty vector if none found
vec<int> findJewelNaive(const Graph &G) {
  vec<int> v;

  while (1) {
    nextPathInPlace(G, v, 4, false, false);

    if (v.empty()) break;

    for (int v5 : G[v[0]]) {
      if (G.areNeighbours(v.back(), v5)) {
        v.push_back(v5);
        if (isJewel(G, v)) return v;
        v.pop_back();
      }
    }
  }

  return vec<int>();

  // vec<int> v(5);
  // do {
  //   if (isJewel(G, v)) return v;

  //   nextTupleInPlace(v, G.n);
  // } while (!isAllZeros(v));

  // return vec<int>();
}

bool containsJewelNaive(const Graph &G) {
  auto v = findJewelNaive(G);
  return !v.empty();
}

// bool containsJewel(const Graph &G) {
//   for (int v2 = 0; v2 < G.n; v2++) {
//     for (int v3 : G[v2]) {
//       for (int v5 = 0; v5 < G.n; v5++) {
//         if (v5 == v2 || v5 == v3) continue;

//       }
//     }
//   }
// }