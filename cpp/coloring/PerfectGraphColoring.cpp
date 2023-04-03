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

Graph generate_graph(std::optional<NetworKit::Graph> graph) {
    std::map<NetworKit::node, int> vertices;
    int nodes = 0;
    graph->forNodes([&](NetworKit::node v) {
        vertices[v] = nodes++;
    });
    std::vector<std::vector<int>> M(nodes);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        M[vertices[u]].push_back(vertices[v]), M[vertices[v]].push_back(vertices[u]);
    });
    return Graph(M);
}

NetworKit::Graph retrieve_graph(const Graph &G) {
  NetworKit::Graph graph(G.n, false, false);
  for (int i = 0; i < G.n; i++) {
    for (int j = i + 1; j < G.n; j++) {
      if (G.areNeighbours(i, j)) {
        graph.addEdge(i, j);
      }
    }
  }
  return graph;
}

int getTheta(const Graph &G, const std::vector<int> &isNodeRemoved) {
  NetworKit::Graph graph(retrieve_graph(G));
  std::vector<int> keep_nodes(isNodeRemoved.size());
  for (int i = 0; i < isNodeRemoved.size(); i++) {
    keep_nodes[i] = !isNodeRemoved[i];
  }
  std::vector<int> indices(graph.numberOfNodes());
  std::partial_sum(keep_nodes.begin(), keep_nodes.end(), indices.begin());
  if (indices.empty() || indices.back() == 0) {
    return 0;
  }

  int m = 0;
  std::vector<int> from, to;
  graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
      if (keep_nodes[u] && keep_nodes[v]) {
          from.push_back(indices[u]), to.push_back(indices[v]), m++;
          assert(G.areNeighbours(u, v));
      }
  });

  if (indices.back() == 0) {
      return 0;
  }
  if (indices.back() == 1) {
      return 1;
  }
  if (indices.back() == 2) {
      return 1 + (m == 0);
  }

  auto e = std::make_tuple(indices.back(), m, from, to);
  static std::map<std::tuple<int, int, std::vector<int>, std::vector<int>>, int> CACHE;
  if (CACHE.count(e) > 0) {
    return CACHE[e];
  }

  double theta = get_theta(indices.back(), m, from.data(), to.data());
  int theta_int = theta + 0.5;
  double eps = 0.3;
  if (abs(theta - theta_int) > eps) {
    throw std::logic_error("Theta returned non-integer for a Perfect Graph: " + std::to_string(theta));
  }
  return CACHE[e] = theta_int;
}

int getOmega(const Graph &G) {
  return getTheta(G.getComplement(), std::vector<int>(G.n, 0));
}

std::vector<int> getMaxCardStableSet(const Graph &G) {
  int thetaG = getTheta(G, std::vector<int>(G.n, 0));

  std::vector<int> isNodeRemoved(G.n, 0);
  std::vector<int> res;

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
  std::partial_sum(c.begin(), c.end(), c.begin());
  std::vector<std::vector<int>> nneighbors(c.back());
  for (int i = 0; i < G.n; i++) {
    for (int j = i + 1; j < G.n; j++) {
      if (G.areNeighbours(i, j)) {
        for (int ni = i == 0 ? 0 : c[i - 1]; ni < c[i]; ni++) {
          for (int nj = c[j - 1]; nj < c[j]; nj++) {
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
    while (wsk < G.n && c[wsk] <= nSSnode) wsk++;

    if (ret.empty() || ret.back() != wsk) {
      ret.push_back(wsk);
    }
  }
  return ret;
}

namespace Koala {

int PerfectGraphColoring::get_omega() {
    return getOmega(generate_graph(graph));
}

int PerfectGraphColoring::get_chi() {
    run();
    int out = 0;
    for (auto &[_, v] : colors) {
        out = std::max(out, v);
    }
    return out;
}

std::set<NetworKit::node> PerfectGraphColoring::get_stable_set_intersecting_all_maximum_cliques() {
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
            K.push_back(getMaxCardClique(Gprim));
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
