#include "color.h"
#include <tuple>
#include <map>
#include "commons.h"
// #include "perfect.h"
// #include "testCommons.h"
using namespace std;

#include "theta.h"

std::tuple<int, int, vec<int>, vec<int>> getGraphEdges(const Graph &G, const vec<int> &isNodeRemoved) {
  vec<int> nodeStays(isNodeRemoved.size());
  for (int i = 0; i < isNodeRemoved.size(); i++) {
    nodeStays[i] = !isNodeRemoved[i];
  }

  while (nodeStays.size() < G.n) {
    nodeStays.push_back(1);
  }

  vec<int> nodeNr = getPrefSum(nodeStays);

  vec<int> from;
  vec<int> to;

  if (nodeNr.empty() || nodeNr.back() == 0) {
    return {0, 0, from, to};
  }

  int m = 0;
  for (int i = 0; i < G.n; i++) {
    if (!nodeStays[i]) continue;

    for (auto j : G[i]) {
      if (!nodeStays[j]) continue;

      if (j > i) {
        from.push_back(nodeNr[i]);
        to.push_back(nodeNr[j]);
        m++;
      }
    }
  }

  return {nodeNr.back(), m, from, to};
}

int getTheta(const Graph &G, const vec<int> &isNodeRemoved, bool gatherStats) {
  auto e = getGraphEdges(G, isNodeRemoved);

  if (get<0>(e) == 0) {
    return 0;
  }

  if (get<0>(e) == 1) {
    return 1;
  }

  if (get<0>(e) == 2) {
    if (get<1>(e) == 0)
      return 2;
    else
      return 1;
  }

  static std::map<std::tuple<int, int, vec<int>, vec<int>>, int> MEM;
  if (MEM.count(e) > 0) {
    return MEM[e];
  }

  double th = theta(get<0>(e), get<1>(e), get<2>(e).data(), get<3>(e).data());

  if (th == -1) {
    throw std::logic_error("Theta returned -1");
  }

  int thInt = th + 0.5;

  double eps = 0.3;
  if (abs(th - thInt) > eps) {
    throw std::logic_error("Theta returned non-integer for a Perfect Graph: " + to_string(th));
  }

  MEM[e] = thInt;

  return thInt;
}

int getOmega(const Graph &G, bool gatherStats) {
  return getTheta(G.getComplement(), vec<int>(), gatherStats);
}

bool isStableSet(const Graph &G, vec<int> nodes) {
  if (!isDistinctValues(nodes)) return false;

  for (int i = 0; i < nodes.size(); i++) {
    for (int j = i + 1; j < nodes.size(); j++) {
      if (G.areNeighbours(nodes[i], nodes[j])) return false;
    }
  }

  return true;
}

vec<int> getMaxCardStableSet(const Graph &G, bool gatherStats) {
  int thetaG = getTheta(G, vec<int>(), gatherStats);

  vec<int> isNodeRemoved(G.n, 0);
  vec<int> res;

  for (int i = 0; i < G.n; i++) {
    isNodeRemoved[i] = 1;
    int newTheta = getTheta(G, isNodeRemoved, gatherStats);

    if (newTheta != thetaG) {
      isNodeRemoved[i] = 0;
      res.push_back(i);
    }
  }

  return res;
}

vec<int> getMaxCardClique(const Graph &G, bool gatherStats) {
  return getMaxCardStableSet(G.getComplement(), gatherStats);
}

vec<int> getSSIntersectingCliques(const Graph &G, vec<vec<int>> K, bool gatherStats) {
  vec<int> c(G.n);

  for (auto k : K) {
    for (auto v : k) {
      c[v]++;
    }
  }

  vec<int> prefC = getPrefSum(c);

  vec<vec<int>> nneighbors(prefC.back());

  for (int i = 0; i < G.n; i++) {
    for (int j = i + 1; j < G.n; j++) {
      if (G.areNeighbours(i, j)) {
        for (int ni = i == 0 ? 0 : prefC[i - 1]; ni < prefC[i]; ni++) {
          for (int nj = prefC[j - 1]; nj < prefC[j]; nj++) {
            nneighbors[ni].push_back(nj);
            nneighbors[nj].push_back(ni);
          }
        }
      }
    }
  }

  Graph nG(nneighbors);

  vec<int> nSS = getMaxCardStableSet(Graph(nneighbors), gatherStats);

  vec<int> ret;
  int wsk = 0;
  for (int nSSnode : nSS) {
    while (wsk < G.n && prefC[wsk] <= nSSnode) wsk++;

    if (ret.empty() || ret.back() != wsk) {
      ret.push_back(wsk);
    }
  }

  return ret;
}

vec<int> getSSIntersectingAllMaxCardCliques(const Graph &G, bool gatherStats) {
  vec<vec<int>> K;
  K.push_back(getMaxCardClique(G, gatherStats));

  int omegaG = getOmega(G, gatherStats);

  while (true) {
    vec<int> S = getSSIntersectingCliques(G, K, gatherStats);

    vec<int> compS = getComplementNodesVec(G.n, S);

    Graph Gprim = G.getInducedStrong(compS);

    if (getOmega(Gprim, gatherStats) < omegaG) {
      return S;
    } else {
      K.push_back(getMaxCardClique(Gprim, gatherStats));

      for (int i = 0; i < K.back().size(); i++) {
        K.back()[i] = compS[K.back()[i]];
      }
    }
  }
}

vec<int> color(const Graph &G, bool gatherStats) {
  if (G.n == 0) return vec<int>();
  if (G.n == 1) return vec<int>{0};

  vec<int> ret(G.n);

  vec<int> SS = getSSIntersectingAllMaxCardCliques(G, gatherStats);
  vec<int> compSS = getComplementNodesVec(G.n, SS);

  vec<int> colorCompSS = color(G.getInducedStrong(compSS));
  for (int i = 0; i < compSS.size(); i++) {
    ret[compSS[i]] = colorCompSS[i] + 1;
  }

  return ret;
}