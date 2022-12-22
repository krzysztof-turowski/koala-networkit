/*
 * NearCleaners.cpp
 *
 *  Created on: 22.12.2022
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <cassert>
#include <set>

#include <boost/dynamic_bitset.hpp>

#include <traversal/PathInplace.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

namespace Koala {

static const NetworKit::count MAX_LENGTH = 512;
using Bitset = boost::dynamic_bitset<uint64_t>;

auto get_all_paths(const NetworKit::Graph &graph, NetworKit::count length) {
  std::vector<std::vector<NetworKit::node>> out;
  std::vector<NetworKit::node> P;
  while (Koala::Traversal::NextPathInplace(
          graph, length, P, Koala::Traversal::PathInplaceMode::INDUCED_PATH)) {
      out.push_back(P);
  }
  return out;
}

template <typename Container>
Bitset get_bitset(NetworKit::count length, Container positions) {
    if (length >= MAX_LENGTH) {
        throw std::length_error("Algorithm cannot be run for graphs on more than 512 vertices");
    }
    Bitset out(length);
    for (auto i : positions) {
        out.set(i);
    }
    return out;
}

auto all_shortest_paths_with_penultimate(const NetworKit::Graph &graph, auto test) {
  unsigned n = graph.upperNodeIdBound(), infinity = std::numeric_limits<unsigned>::max();
  std::vector<std::vector<NetworKit::count>> D(n, std::vector<NetworKit::count>(n, infinity));
  std::vector<std::vector<NetworKit::node>> penultimate(
      n, std::vector<NetworKit::node>(n, NetworKit::none));

  graph.forNodes([&](NetworKit::node i) { D[i][i] = 0; });
  graph.forEdges([&](NetworKit::node i, NetworKit::node j) {
      D[i][j] = D[j][i] = 1, penultimate[i][j] = i, penultimate[j][i] = j;
  });
  for (auto k : graph.nodeRange()) {
      if (!test(k)) {
          continue;
      }
      for (auto i : graph.nodeRange()) {
          if (i == k) {
              continue;
          }
          for (auto j : graph.nodeRange()) {
              if (j != i && j != k && D[i][j] > D[i][k] + D[k][j]) {
                D[i][j] = D[i][k] + D[k][j], penultimate[i][j] = penultimate[k][j];
              }
          }
      }
  }
  return std::make_tuple(D, penultimate);
}

bool check_odd_hole_with_near_cleaner(
          const NetworKit::Graph &graph, const Bitset &S,
          const std::vector<std::vector<NetworKit::node>> &triplePaths) {
    auto infinity = std::numeric_limits<NetworKit::count>::max();
    auto [D, penultimate] = all_shortest_paths_with_penultimate(
        graph, [&](auto v) { return !S.test(v); });
    for (const auto &y1 : graph.nodeRange()) {
        if (S.test(y1)) {
            continue;
        }
        for (const auto &triple : triplePaths) {
            if (std::ranges::any_of(
                    std::views::iota(0, 3), [&](auto i) { return triple[i] == y1; })) {
                continue;
            }
            auto x1 = triple[0], x3 = triple[1], x2 = triple[2];
            if (D[x1][y1] == infinity || D[x2][y1] == infinity) {
                continue;
            }
            auto y2 = penultimate[x2][y1];
            auto n = D[x2][y1];
            if (D[x1][y1] + 1 != n || D[x1][y2] != n || D[x3][y1] < n || D[x3][y2] < n) {
                continue;
            }
            return true;
        }
    }
    return false;
}

inline bool is_relevant_triple(
        const NetworKit::Graph &graph, NetworKit::node a, NetworKit::node b, NetworKit::node c) {
    return !(a == b || graph.hasEdge(a, b) || (graph.hasEdge(a, c) && graph.hasEdge(b, c)));
}

Bitset get_x_for_relevant_triple(
        const NetworKit::Graph &graph, NetworKit::node a, NetworKit::node b, NetworKit::node c) {
    auto anticomponents_Nab = PerfectGraphRecognition::getAuxiliaryComponents(graph, {a, b});
    auto non_edge_c = [&](auto v) { return !graph.hasEdge(c, v); };
    unsigned threshold = 0;
    for (const auto &component : anticomponents_Nab) {
        if (component.size() <= threshold) {
            continue;
        }
        bool contains_nonneighbour_of_c = std::any_of(
            component.cbegin(), component.cend(), non_edge_c);
        if (contains_nonneighbour_of_c) {
            threshold = component.size();
        }
    }
    std::vector<NetworKit::node> Y;
    for (const auto &component : anticomponents_Nab) {
        if (component.size() > threshold) {
            Y.insert(Y.end(), component.begin(), component.end());
        }
    }
    auto it = std::find_if(
        anticomponents_Nab.cbegin(), anticomponents_Nab.cend(),
        [&](auto component) {
            return std::any_of(component.cbegin(), component.cend(), non_edge_c);
    });
    std::vector<NetworKit::node> W =
        it == anticomponents_Nab.cend() ? std::vector<NetworKit::node>() : *it;
    W.push_back(c);
    W.insert(W.end(), Y.begin(), Y.end());
    auto Z = PerfectGraphRecognition::getAllCompleteVertices(graph, W);
    return get_bitset(graph.upperNodeIdBound(), Y) | get_bitset(graph.upperNodeIdBound(), Z);
}

auto get_possible_near_cleaners(const NetworKit::Graph &graph) {
    std::vector<Bitset> Ns, Xs;
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        auto C = PerfectGraphRecognition::getAllCompleteVertices(graph, {u, v});
        Ns.push_back(get_bitset(graph.upperNodeIdBound(), C));
    });
    graph.forNodePairs([&](NetworKit::node a, NetworKit::node b) {
        if (!graph.hasEdge(a, b)) {
            for (const auto &c : graph.nodeRange()) {
                if (is_relevant_triple(graph, a, b, c)) {
                    Xs.push_back(get_x_for_relevant_triple(graph, a, b, c));
                }
            }
        }
    });
    std::set<Bitset> out;
    for (const auto &N : Ns) {
        for (const auto &X : Xs) {
            out.insert(X | N);
        }
    }
    return out;
}

bool PerfectGraphRecognition::containsNearCleanerOddHole(const NetworKit::Graph &graph) {
    auto triplePaths = get_all_paths(graph, 3);
    for (const auto &X : get_possible_near_cleaners(graph)) {
        if (check_odd_hole_with_near_cleaner(graph, X, triplePaths)) {
            return true;
        }
    }
    return false;
}

} /* namespace Koala */
