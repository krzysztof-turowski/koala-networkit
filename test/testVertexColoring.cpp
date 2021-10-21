#include <gtest/gtest.h>

#include <list>

#include <coloring/GreedyVertexColoring.hpp>

struct VertexColoringParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int colors;
};

class GreedyVertexColoringRandomSequentialTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class GreedyVertexColoringLargestFirstTest
    : public testing::TestWithParam<VertexColoringParameters> { };

NetworKit::Graph build_graph(const int &N, const std::list<std::pair<int, int>> &E) {
  NetworKit::Graph G(N, false, false);
  for (const auto &[u, v] : E) {
    G.addEdge(u, v);
  }
  return G;
}

TEST_P(GreedyVertexColoringRandomSequentialTest, test) {
  VertexColoringParameters const& parameters = GetParam();
  NetworKit::Graph G = build_graph(parameters.N, parameters.E);
  std::map<NetworKit::node, int> S;
  EXPECT_EQ(Koala::GreedyVertexColoring().randomSequential(G, S), parameters.colors);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, GreedyVertexColoringRandomSequentialTest, testing::Values(
        VertexColoringParameters{4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, 2}
));

TEST_P(GreedyVertexColoringLargestFirstTest, test) {
  VertexColoringParameters const& parameters = GetParam();
  NetworKit::Graph G = build_graph(parameters.N, parameters.E);
  std::map<NetworKit::node, int> S;
  EXPECT_EQ(Koala::GreedyVertexColoring().largestFirst(G, S), parameters.colors);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, GreedyVertexColoringLargestFirstTest, testing::Values(
        VertexColoringParameters{
            7, {{0, 2}, {0, 3}, {0, 5}, {0, 6}, {1, 2}, {1, 4}, {1, 5}, {1, 6}, {2, 3}, {2, 4},
                {3, 4}}, 4},
        VertexColoringParameters{
            8, {{0, 1}, {0, 2}, {0, 4}, {0, 6}, {0, 7}, {1, 2}, {1, 3}, {1, 7}, {2, 3}, {2, 4},
                {3, 5}, {3, 7}, {4, 5}, {4, 6}, {5, 6}, {5, 7}, {6, 7}}, 5},
        VertexColoringParameters{
            10, {{0, 2}, {0, 3}, {0, 4}, {0, 5}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 6}, {2, 8},
                 {3, 4}, {3, 6}, {4, 6}, {5, 6}, {5, 7}, {5, 9}, {7, 8}, {7, 9}, {8, 9}}, 4},
        VertexColoringParameters{
            10, {{0, 1}, {0, 2}, {0, 3}, {0, 5}, {0, 6}, {0, 7}, {1, 2}, {1, 3}, {1, 4}, {1, 6},
                 {1, 7}, {2, 3}, {2, 4}, {2, 5}, {2, 7}, {3, 4}, {3, 5}, {3, 6}, {4, 5}, {4, 7},
                 {4, 8}, {4, 9}, {5, 6}, {5, 8}, {5, 9}, {6, 7}, {6, 8}, {6, 9}, {7, 8}, {7, 9},
                 {8, 9}}, 6},
        VertexColoringParameters{
            10, {{0, 1}, {0, 4}, {0, 5}, {0, 6}, {0, 8}, {1, 2}, {1, 4}, {1, 5}, {1, 7}, {1, 8},
                 {2, 3}, {2, 4}, {2, 6}, {2, 7}, {2, 9}, {3, 5}, {3, 6}, {3, 7}, {3, 8}, {4, 8},
                 {5, 6}, {5, 7}, {5, 8}, {5, 9}, {6, 9}, {7, 9}, {8, 9}}, 5}
));
