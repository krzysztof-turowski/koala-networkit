#include <algorithm>
#include <limits>
#include <vector>

#include <structures/dynamic_tree/NaiveDynamicTree.hpp>

namespace Koala {

std::vector<NetworKit::node> getPath(
    const NetworKit::Graph &graph, NetworKit::node u, NetworKit::node v);

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
  std::vector<NetworKit::node> s;
  std::vector<bool> visited(graph.numberOfNodes(), false);
  s.push_back(v);
  while (!s.empty()) {
    NetworKit::node x = s.back();
    s.pop_back();
    visited[x] = true;
    root = std::min(root, x);
    graph.forNeighborsOf(x, [&](NetworKit::node y) {
      if (!visited[y]) {
        s.push_back(y);
      }
    });
  }
  return root;
}

template <typename Value>
void NaiveDynamicTree<Value>::pathAdd(NetworKit::node u, NetworKit::node v, Value c) {
  auto path = getPath(graph, u, v);
  for (std::size_t i = 1; i < path.size(); ++i) {
    weights[path[i - 1]][path[i]] -= c;
    weights[path[i]][path[i - 1]] += c;
  }
}

template <typename Value>
NetworKit::Edge NaiveDynamicTree<Value>::pathMin(NetworKit::node u, NetworKit::node v) {
  auto path = getPath(graph, u, v);
  Value m = std::numeric_limits<Value>::max();
  NetworKit::Edge minimumEdge;
  for (std::size_t i = 1; i < path.size(); ++i) {
    Value w = weights[path[i - 1]][path[i]];
    if (w >= 0 && w <= m) {
      m = w;
      minimumEdge = NetworKit::Edge(path[i - 1], path[i]);
    } else if (w < 0 && Value{1} + w <= m) {
      m = Value{1} + w;
      minimumEdge = NetworKit::Edge(path[i - 1], path[i]);
    }
  }
  return minimumEdge;
}

template <typename Value>
Value NaiveDynamicTree<Value>::pathSum(NetworKit::node u, NetworKit::node v) {
  auto path = getPath(graph, u, v);
  Value s = 0;
  for (std::size_t i = 1; i < path.size(); ++i) {
    s += weights[path[i - 1]][path[i]];
  }
  return s;
}

template <typename Value>
Value NaiveDynamicTree<Value>::getWeight(NetworKit::node u, NetworKit::node v) const {
  return weights[u][v];
}

template <typename Value>
void NaiveDynamicTree<Value>::addWeight(NetworKit::node u, NetworKit::node v, Value value) {
  weights[u][v] += value;
}

std::vector<NetworKit::node> getPath(
    const NetworKit::Graph &graph, NetworKit::node u, NetworKit::node v) {
  std::vector<NetworKit::node> s;
  std::vector<NetworKit::node> parent(graph.numberOfNodes(), NetworKit::none);
  s.push_back(u);
  parent[u] = u;
  while (!s.empty()) {
    NetworKit::node x = s.back();
    if (x == v) {
      break;
    }
    s.pop_back();
    graph.forNeighborsOf(x, [&](NetworKit::node y) {
      if (parent[y] == NetworKit::none) {
        parent[y] = x;
        s.push_back(y);
      }
    });
  }
  std::vector<NetworKit::node> path;
  NetworKit::node x = v;
  while (x != u) {
    path.push_back(x);
    x = parent[x];
  }
  path.push_back(u);

  std::reverse(path.begin(), path.end());

  return path;
}

template class NaiveDynamicTree<int>;
template class NaiveDynamicTree<double>;

}  // namespace Koala
