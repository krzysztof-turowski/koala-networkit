#include "commons.h"
#include <map>
#include <queue>
#include <set>
#include <unordered_set>

void Graph::calculateNeighboursLists() {
  _neighbours.clear();
  for (int i = 0; i < n; i++) {
    _neighbours.push_back(std::vector<int>());
    for (int j = 0; j < n; j++) {
      if (areNeighbours(i, j)) _neighbours[i].push_back(j);
    }
  }
}

void Graph::calculateFirstNextNeighbours() {
  _first_neighbour = std::vector<int>(n, -1);
  _next_neighbour = std::vector<std::vector<int>>(n);
  for (int i = 0; i < n; i++) {
    _next_neighbour[i] = std::vector<int>(n, -2);
  }

  for (int i = 0; i < n; i++) {
    if (!_neighbours[i].empty()) {
      _first_neighbour[i] = _neighbours[i].front();
      _next_neighbour[i][_neighbours[i].back()] = -1;
    }

    for (int j = 1; j < _neighbours[i].size(); j++) {
      _next_neighbour[i][_neighbours[i][j - 1]] = _neighbours[i][j];
    }
  }
}

void Graph::checkSymmetry() {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (_matrix[i][j] != _matrix[j][i]) {
        throw std::invalid_argument("Graph is not symmetrical. (" + std::to_string(i) + ", " + std::to_string(j) + ").");
      }
    }
  }
}

Graph::Graph(int n) : n(n), _neighbours(n), _matrix(n) {
  for (int i = 0; i < n; i++) {
    _matrix[i].resize(n);
  }

  calculateFirstNextNeighbours();
}

Graph::Graph(std::vector<std::vector<int>> neighbours) : n(neighbours.size()), _neighbours(neighbours) {
  _matrix.resize(n);
  for (int i = 0; i < n; i++) {
    _matrix[i].resize(n);
  }

  for (int i = 0; i < n; i++) {
    for (int j : _neighbours[i]) {
      if (j < 0 || j >= n) {
        throw std::invalid_argument("Graph initialization from neighbours array failed. Neighbour out of range.");
      }
      _matrix[i][j] = true;
    }
  }

  calculateFirstNextNeighbours();
  checkSymmetry();
}

int Graph::getFirstNeighbour(int a) const { return _first_neighbour[a]; }

int Graph::getNextNeighbour(int a, int b) const {
  int ret = _next_neighbour[a][b];

  if (ret == -2) {
    char buff[100];
    snprintf(buff, sizeof(buff), "Graph getNextNeighbour failed. %d is not a neighbour of %d.", b, a);

    throw std::invalid_argument(buff);
  }

  return ret;
}

Graph Graph::getComplement() const {
  Graph ret(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) ret._matrix[i][j] = !_matrix[i][j];
    }
  }

  ret.calculateNeighboursLists();
  ret.calculateFirstNextNeighbours();
  return ret;
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

  ret.calculateNeighboursLists();
  ret.calculateFirstNextNeighbours();
  return ret;
}

std::vector<int> getPrefSum(const std::vector<int> &v) {
  if (v.empty()) return std::vector<int>();

  std::vector<int> ret(v.size());
  ret[0] = v[0];
  for (int i = 1; i < v.size(); i++) {
    ret[i] = ret[i - 1] + v[i];
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
