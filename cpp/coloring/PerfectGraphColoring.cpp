/*
 * PerfectGraphColoring.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <coloring/PerfectGraphColoring.hpp>

#include <map>
#include <optional>
#include <tuple>

#include "perfect/commons.h"
#include "perfect/theta.h"

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

int getTheta(const Graph &G, const vec<int> &isNodeRemoved) {
  auto e = getGraphEdges(G, isNodeRemoved);

  if (std::get<0>(e) == 0) {
    return 0;
  }

  if (std::get<0>(e) == 1) {
    return 1;
  }

  if (std::get<0>(e) == 2) {
    if (std::get<1>(e) == 0)
      return 2;
    else
      return 1;
  }

  static std::map<std::tuple<int, int, vec<int>, vec<int>>, int> MEM;
  if (MEM.count(e) > 0) {
    return MEM[e];
  }

  // std::cout << "nodes " << get<0>(e) << std::endl;
  // std::cout << "m " << get<1>(e) << std::endl;
  // std::cout << "isNodeRemoved " << isNodeRemoved.size() << std::endl;
  // for (int i = 0; i < from.size(); i++) {
    // std::cout << "edge " << i << " " << from[i] << " " << to[i] << std::endl;
  // }
  double th = theta(std::get<0>(e), std::get<1>(e), std::get<2>(e).data(), std::get<3>(e).data());

  if (th == -1) {
    throw std::logic_error("Theta returned -1");
  }

  int thInt = th + 0.5;

  double eps = 0.3;
  if (abs(th - thInt) > eps) {
    throw std::logic_error("Theta returned non-integer for a Perfect Graph: " + std::to_string(th));
  }

  MEM[e] = thInt;

  return thInt;
}

int getOmega(const Graph &G) {
  return getTheta(G.getComplement(), vec<int>());
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

vec<int> getMaxCardStableSet(const Graph &G) {
  int thetaG = getTheta(G, vec<int>());

  vec<int> isNodeRemoved(G.n, 0);
  vec<int> res;

  for (int i = 0; i < G.n; i++) {
    isNodeRemoved[i] = 1;
    int newTheta = getTheta(G, isNodeRemoved);

    if (newTheta != thetaG) {
      isNodeRemoved[i] = 0;
      res.push_back(i);
    }
  }

  return res;
}

vec<int> getMaxCardClique(const Graph &G) {
  return getMaxCardStableSet(G.getComplement());
}

vec<int> getSSIntersectingCliques(const Graph &G, vec<vec<int>> K) {
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

  vec<int> nSS = getMaxCardStableSet(Graph(nneighbors));

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

vec<int> getSSIntersectingAllMaxCardCliques(const Graph &G) {
  vec<vec<int>> K;
  K.push_back(getMaxCardClique(G));

  int omegaG = getOmega(G);

  while (true) {
    vec<int> S = getSSIntersectingCliques(G, K);

    vec<int> compS = getComplementNodesVec(G.n, S);

    Graph Gprim = G.getInducedStrong(compS);

    if (getOmega(Gprim) < omegaG) {
      return S;
    } else {
      K.push_back(getMaxCardClique(Gprim));

      for (int i = 0; i < K.back().size(); i++) {
        K.back()[i] = compS[K.back()[i]];
      }
    }
  }
}

vec<int> color(const Graph &G) {
  if (G.n == 0) return vec<int>();
  if (G.n == 1) return vec<int>{0};

  vec<int> ret(G.n);

  vec<int> SS = getSSIntersectingAllMaxCardCliques(G);
  vec<int> compSS = getComplementNodesVec(G.n, SS);

  vec<int> colorCompSS = color(G.getInducedStrong(compSS));
  for (int i = 0; i < compSS.size(); i++) {
    ret[compSS[i]] = colorCompSS[i] + 1;
  }

  return ret;
}

namespace Koala {

void PerfectGraphColoring::run() {
    vec<vec<int>> M;
    M.resize(graph->numberOfNodes());
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        M[u].push_back(v), M[v].push_back(u);
        std::cout << u << " " << v << std::endl;
    });
    Graph G(M);
    auto C = color(G);
    hasRun = true;
}

} /* namespace Koala */
