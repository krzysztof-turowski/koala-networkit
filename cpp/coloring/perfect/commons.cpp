#include "commons.h"
#include <map>
#include <queue>
#include <set>
#include <unordered_set>

Graph::Graph(int n) : n(n), _matrix(n) {
  for (int i = 0; i < n; i++) {
    _matrix[i].resize(n);
  }
}

Graph::Graph(std::vector<std::vector<int>> neighbours) : n(neighbours.size()) {
  _matrix.resize(n);
  for (int i = 0; i < n; i++) {
    _matrix[i].resize(n);
  }

  for (int i = 0; i < n; i++) {
    for (int j : neighbours[i]) {
      if (j < 0 || j >= n) {
        throw std::invalid_argument("Graph initialization from neighbours array failed. Neighbour out of range.");
      }
      _matrix[i][j] = true;
    }
  }
}

Graph Graph::getInducedStrong(std::vector<int> X) const {
  Graph ret(X.size());
  for (int i = 0; i < X.size(); i++) {
    for (int j = 0; j < X.size(); j++) {
      if (areNeighbours(X[i], X[j])) {
        ret._matrix[i][j] = true;
      }
    }
  }
  return ret;
}

std::vector<int> getComplementNodesVec(int n, const std::vector<int> &X) {
  for (int i = 1; i < X.size(); i++) {
    if (X[i - 1] >= X[i]) {
      throw std::invalid_argument("X for getComplementNodesVec should be sorted");
    }
  }

  std::vector<int> res;

  res.reserve(n - X.size());
  int wsk = 0;
  for (int i = 0; i < n; i++) {
    while (wsk < X.size() && X[wsk] < i) wsk++;

    if (wsk >= X.size() || X[wsk] != i) res.push_back(i);
  }

  return res;
}
