/*
 * PerfectGraphColoring.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <map>
#include <optional>
#include <tuple>

#include <coloring/PerfectGraphColoring.hpp>

#include "perfect/commons.h"
#include "perfect/theta.h"

int getTheta(const Graph &G, const std::set<int> &isNodeRemoved) {
    if (G.n - isNodeRemoved.size() == 0) {
        return 0;
    }
    if (G.n - isNodeRemoved.size() == 1) {
        return 1;
    }
    if (G.n - isNodeRemoved.size() == 2) {
        return 1 + G.areNeighbours(0, 1);
    }

  std::map<int, int> V;
  int n = 1, m = 0;
  for (int i = 0; i < G.n; i++) {
      if (!isNodeRemoved.count(i)) {
          V[i] = n++;
      }
  }
  std::vector<int> from, to;
  for (int i = 0; i < G.n; i++) {
      if (isNodeRemoved.count(i)) {
          continue;
      }
      for (auto j : G[i]) {
          if (isNodeRemoved.count(i)) {
              continue;
          }
          if (j > i) {
            from.push_back(V[i]);
            to.push_back(V[j]);
            m++;
          }
      }
  }
  double th = theta(G.n - isNodeRemoved.size(), m, from.data(), to.data());

  if (th == -1) {
    throw std::logic_error("Theta returned -1");
  }
  if (abs(th - static_cast<int>(th + 0.5)) > 0.3) {
    throw std::logic_error("Theta returned non-integer for a Perfect Graph: " + std::to_string(th));
  }
  return static_cast<int>(th + 0.5);
}

int getOmega(const Graph &G) {
  return getTheta(G.getComplement(), std::set<int>());
}

std::vector<int> getMaxCardStableSet(const Graph &G) {
  int thetaG = getTheta(G, std::set<int>());

  std::set<int> isNodeRemoved;
  std::vector<int> res;

  for (int i = 0; i < G.n; i++) {
    isNodeRemoved.insert(i);
    int newTheta = getTheta(G, isNodeRemoved);

    if (newTheta != thetaG) {
      isNodeRemoved.erase(i);
      res.push_back(i);
    }
  }

  return res;
}

std::vector<int> getMaxCardClique(const Graph &G) {
  return getMaxCardStableSet(G.getComplement());
}

std::vector<int> getSSIntersectingCliques(const Graph &G, std::vector<std::vector<int>> K) {
  std::vector<int> c(G.n);

  for (auto k : K) {
    for (auto v : k) {
      c[v]++;
    }
  }

  std::vector<int> prefC = getPrefSum(c);

  std::vector<std::vector<int>> nneighbors(prefC.back());

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

  std::vector<int> nSS = getMaxCardStableSet(Graph(nneighbors));

  std::vector<int> ret;
  int wsk = 0;
  for (int nSSnode : nSS) {
    while (wsk < G.n && prefC[wsk] <= nSSnode) wsk++;

    if (ret.empty() || ret.back() != wsk) {
      ret.push_back(wsk);
    }
  }

  return ret;
}

Graph generate_graph(std::optional<NetworKit::Graph> graph) {
    std::map<NetworKit::node, int> vertices;
    int nodes = 0;
    graph->forNodes([&](NetworKit::node v) {
        vertices[v] = nodes++;
    });
    std::vector<std::vector<int>> M;
    M.resize(nodes);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        M[vertices[u]].push_back(vertices[v]), M[vertices[v]].push_back(vertices[u]);
        std::cout << "EDGE " << u << " " << v << std::endl;
    });
    return Graph(M);
}

namespace Koala {

std::set<NetworKit::node> PerfectGraphColoring::get_stable_set_intersecting_all_maximum_cliques() {
    std::cout << "GRAPH SIZE: " << graph->numberOfNodes() << std::endl;
    Graph G(generate_graph(graph));
    std::vector<int> vertices;
    graph->forNodes([&](NetworKit::node v) {
        vertices.push_back(v);
    });
    std::vector<std::vector<int>> K;
    K.push_back(getMaxCardClique(G));
    int omegaG = K.back().size();
    while (true) {
        std::vector<int> S = getSSIntersectingCliques(G, K);
        std::vector<int> compS = getComplementNodesVec(G.n, S);
        Graph Gprim = G.getInducedStrong(compS);
        if (getOmega(Gprim) < omegaG) {
            std::set<NetworKit::node> out;
            for (auto i : S) {
                out.insert(vertices[i]);
            }
            return out;
        } else {
            K.push_back();

            for (int i = 0; i < K.back().size(); i++) {
              K.back()[i] = compS[K.back()[i]];
            }
        }
    }
}

void PerfectGraphColoring::run() {
    for (int color = 1; graph->numberOfNodes() > 0; color++) {
        for (const auto &v : get_stable_set_intersecting_all_maximum_cliques()) {
            colors[v] = color;
            graph->removeNode(v);
        }
    }
    hasRun = true;
}

} /* namespace Koala */
