#pragma once

#include <algorithm>
#include <limits>
#include <vector>

#include <structures/DynamicTree.hpp>

#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

template <typename Value>
class NaiveDynamicTree final : public DynamicTree<Value> {
 public:
  NaiveDynamicTree(NetworKit::count n, std::vector<std::vector<Value>> &weights);

  void link(NetworKit::node u, NetworKit::node v, Value value) override;
  void cut(NetworKit::node u, NetworKit::node v) override;
  NetworKit::node findRoot(NetworKit::node v) override;
  void pathAdd(NetworKit::node u, NetworKit::node v, Value value);
  NetworKit::Edge pathMin(NetworKit::node u, NetworKit::node v);
  Value pathSum(NetworKit::node u, NetworKit::node v);

  Value getWeight(NetworKit::node u, NetworKit::node v) const;
  void addWeight(NetworKit::node u, NetworKit::node v, Value value);

 private:
  NetworKit::Graph graph;
  std::vector<std::vector<Value>> &weights;
};

template <typename Value>
NaiveDynamicTree<Value>::NaiveDynamicTree(
    NetworKit::count n, std::vector<std::vector<Value>> &weights)
    : graph(NetworKit::Graph(n)), weights(weights) {}

template <typename Value>
void NaiveDynamicTree<Value>::link(NetworKit::node u, NetworKit::node v, Value) {
  graph.addEdge(u, v);
}

template <typename Value>
void NaiveDynamicTree<Value>::cut(NetworKit::node u, NetworKit::node v) {
  if (graph.hasEdge(u, v)) {
    graph.removeEdge(u, v);
  }
}

template <typename Value>
NetworKit::node NaiveDynamicTree<Value>::findRoot(NetworKit::node v) {
  NetworKit::node root = v;
  std::vector<NetworKit::node> stack{v};
  std::vector<bool> visited(graph.numberOfNodes(), false);
  while (!stack.empty()) {
    NetworKit::node u = stack.back();
    stack.pop_back();
    visited[u] = true;
    root = std::min(root, u);
    graph.forNeighborsOf(u, [&](NetworKit::node neighbor) {
      if (!visited[neighbor]) {
        stack.push_back(neighbor);
      }
    });
  }
  return root;
}

inline std::vector<NetworKit::node> getPath(
    const NetworKit::Graph &graph, NetworKit::node u, NetworKit::node v) {
  std::vector<NetworKit::node> stack{u};
  std::vector<NetworKit::node> parent(graph.numberOfNodes(), NetworKit::none);
  parent[u] = u;
  while (!stack.empty()) {
    NetworKit::node current = stack.back();
    if (current == v) {
      break;
    }
    stack.pop_back();
    graph.forNeighborsOf(current, [&](NetworKit::node neighbor) {
      if (parent[neighbor] == NetworKit::none) {
        parent[neighbor] = current;
        stack.push_back(neighbor);
      }
    });
  }
  std::vector<NetworKit::node> path;
  for (NetworKit::node current = v; current != u; current = parent[current]) {
    path.push_back(current);
  }
  path.push_back(u);
  std::reverse(path.begin(), path.end());
  return path;
}

template <typename Value>
void NaiveDynamicTree<Value>::pathAdd(NetworKit::node u, NetworKit::node v, Value value) {
  auto path = getPath(graph, u, v);
  for (std::size_t i = 1; i < path.size(); ++i) {
    weights[path[i - 1]][path[i]] -= value;
    weights[path[i]][path[i - 1]] += value;
  }
}

template <typename Value>
NetworKit::Edge NaiveDynamicTree<Value>::pathMin(NetworKit::node u, NetworKit::node v) {
  auto path = getPath(graph, u, v);
  Value minimum = std::numeric_limits<Value>::max();
  NetworKit::Edge minimumEdge;
  for (std::size_t i = 1; i < path.size(); ++i) {
    Value weight = weights[path[i - 1]][path[i]];
    if (weight >= 0 && weight <= minimum) {
      minimum = weight;
      minimumEdge = NetworKit::Edge(path[i - 1], path[i]);
    } else if (weight < 0 && Value{1} + weight <= minimum) {
      minimum = Value{1} + weight;
      minimumEdge = NetworKit::Edge(path[i - 1], path[i]);
    }
  }
  return minimumEdge;
}

template <typename Value>
Value NaiveDynamicTree<Value>::pathSum(NetworKit::node u, NetworKit::node v) {
  auto path = getPath(graph, u, v);
  Value sum = 0;
  for (std::size_t i = 1; i < path.size(); ++i) {
    sum += weights[path[i - 1]][path[i]];
  }
  return sum;
}

template <typename Value>
Value NaiveDynamicTree<Value>::getWeight(NetworKit::node u, NetworKit::node v) const {
  return weights[u][v];
}

template <typename Value>
void NaiveDynamicTree<Value>::addWeight(
    NetworKit::node u, NetworKit::node v, Value value) {
  weights[u][v] += value;
}

};  // namespace Koala
