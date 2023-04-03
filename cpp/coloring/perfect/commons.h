#pragma once

#include <algorithm>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

struct Graph {
  int n;
  explicit Graph(int n);
  explicit Graph(std::vector<std::vector<int>> neighbours);

  // List of neighbors is not necessarily sorted
  std::vector<int> &operator[](int index) { return _neighbours[index]; }
  const std::vector<int> &operator[](int index) const { return _neighbours[index]; }
  bool operator==(const Graph &c) const { return _matrix == c._matrix; }
  bool areNeighbours(int a, int b) const { return _matrix[a][b]; }
  // Returns first neighbor of a. Returns -1 if a has no neighbors.
  int getFirstNeighbour(int a) const;
  // Returns next neighbor of a, after b. Returns -1 if b is last, throws invalid_argument if b is not a
  // neighbor of a. Guaranteed to be consistent with G[a] ordering.
  int getNextNeighbour(int a, int b) const;

  Graph getComplement() const;
  // Returns G' - Graph induced on X. Set of vertices is renamed to match X. e.g. if only node 0 is removed,
  // all other nodes will have number lowered by 1.
  Graph getInducedStrong(std::vector<int> X) const;

 private:
  std::vector<std::vector<int>> _neighbours;
  std::vector<std::vector<int>> _matrix;
  std::vector<int> _first_neighbour;
  std::vector<std::vector<int>> _next_neighbour;
  void checkSymmetry();
  void calculateNeighboursLists();
  // Assumes _neighbours has been calculated;
  void calculateFirstNextNeighbours();
};

std::vector<int> getComplementNodesVec(int n, const std::vector<int> &X);
