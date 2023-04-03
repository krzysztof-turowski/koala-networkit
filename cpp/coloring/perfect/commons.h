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

  bool areNeighbours(int a, int b) const { return _matrix[a][b]; }
  // Returns G' - Graph induced on X. Set of vertices is renamed to match X. e.g. if only node 0 is removed,
  // all other nodes will have number lowered by 1.
  Graph getInducedStrong(std::vector<int> X) const;

 private:
  std::vector<std::vector<int>> _matrix;
};

std::vector<int> getComplementNodesVec(int n, const std::vector<int> &X);
