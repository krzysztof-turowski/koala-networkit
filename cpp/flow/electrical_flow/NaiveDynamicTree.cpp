#include <algorithm>
#include <cassert>
#include <stack>
#include <utility>
#include <vector>

#include <flow/electrical_flow/DynamicTree.hpp>

namespace Koala {

std::vector<int> getPath(const NetworKit::Graph &graph, int u, int v);

DynamicTree::DynamicTree(int n, std::vector<std::vector<double>> &weights)
    : graph(NetworKit::Graph(n)), weights(weights) {}

void DynamicTree::link(int u, int v) { graph.addEdge(u, v); }

void DynamicTree::cut(int u, int v) {
  if (graph.hasEdge(u, v)) {
    graph.removeEdge(u, v);
  }
}

int DynamicTree::findRoot(int v) {
  int root = v;
  std::vector<int> s;
  std::vector<bool> visited(graph.numberOfNodes(), false);
  s.push_back(v);
  while (!s.empty()) {
    int x = s.back();
    s.pop_back();
    visited[x] = true;
    root = std::min(root, x);
    graph.forNeighborsOf(x, [&](int y) {
      if (!visited[y]) {
        s.push_back(y);
      }
    });
  }
  return root;
}

void DynamicTree::pathAdd(int u, int v, double c) {
  auto path = getPath(graph, u, v);
  for (int i = 1; i < path.size(); ++i) {
    weights[path[i - 1]][path[i]] -= c;
    weights[path[i]][path[i - 1]] += c;
  }
}

std::pair<int, int> DynamicTree::pathMin(int u, int v) {
  auto path = getPath(graph, u, v);
  double m = 1e+37;
  std::pair<int, int> mv = {-1, -1};
  for (int i = 1; i < path.size(); ++i) {
    double w = weights[path[i - 1]][path[i]];
    if (w >= 0 && w <= m) {
      m = w;
      mv = {path[i - 1], path[i]};
    } else if (w < 0 && 1.0 + w <= m) {
      m = 1.0 + w;
      mv = {path[i - 1], path[i]};
    }
  }
  return mv;
}

double DynamicTree::pathSum(int u, int v) {
  auto path = getPath(graph, u, v);
  double s = 0;
  for (int i = 1; i < path.size(); ++i) {
    s += weights[path[i - 1]][path[i]];
  }
  return s;
}

std::vector<int> getPath(const NetworKit::Graph &graph, int u, int v) {
  std::vector<int> s;
  std::vector<int> parent(graph.numberOfNodes(), -1);
  s.push_back(u);
  parent[u] = u;
  while (!s.empty()) {
    int x = s.back();
    if (x == v) {
      break;
    }
    s.pop_back();
    graph.forNeighborsOf(x, [&](int y) {
      if (parent[y] == -1) {
        parent[y] = x;
        s.push_back(y);
      }
    });
  }
  std::vector<int> path;
  int x = v;
  while (x != u) {
    path.push_back(x);
    x = parent[x];
  }
  path.push_back(u);

  reverse(path.begin(), path.end());

  return path;
}

}  // namespace Koala
