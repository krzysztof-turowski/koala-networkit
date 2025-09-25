#pragma once

#include <set>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

#include <matching/gaussian_matching/utils.hpp>

typedef std::set<std::pair<int, int>> Matching;

namespace Koala {
class BipartiteGaussianMatching : public NetworKit::Algorithm {
 public:
  explicit BipartiteGaussianMatching(const NetworKit::Graph &G);
  void run();
  Matching getMatching();

  NetworKit::Graph G;
  MatZp AG;
  Matching M;

  std::vector<int> U, V;
  std::vector<int> bpIdx;
  std::vector<int> oldIdx;
};
}  // namespace Koala
