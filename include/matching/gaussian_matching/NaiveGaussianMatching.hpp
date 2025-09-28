#pragma once

#include <set>
#include <vector>

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <matching/gaussian_matching/utils.hpp>

typedef std::set<std::pair<int, int>> Matching;

namespace Koala {
class NaiveGaussianMatching : public NetworKit::Algorithm {
 public:
  explicit NaiveGaussianMatching(const NetworKit::Graph &G);
  void run();
  Matching getMatching();

  NetworKit::Graph G;
  Matching M;
};
}  // namespace Koala
