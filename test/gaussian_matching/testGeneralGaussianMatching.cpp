#include <gtest/gtest.h>

#include <list>
#include <vector>

#include "matching/gaussian_matching/GeneralGaussianMatching.hpp"

#include "../helpers.hpp"

class GenTest : public testing::Test {};

TEST(GenTest, testSuccess) {
  const std::list<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {1, 3}};
  auto G = build_graph(4, edges, false);
  int n = G.numberOfNodes();

  Koala::GeneralGaussianMatching gen(G);
  gen.run();

  auto M = gen.getMatching();
  EXPECT_EQ(M.size(), n / 2);

  std::vector<int> counts(n, 0);
  for (auto [u, v] : M) {
    counts[u]++;
    counts[v]++;
    EXPECT_TRUE(G.hasEdge(u, v));
  }

  for (auto c : counts) {
    EXPECT_EQ(c, 1);
  }
}
