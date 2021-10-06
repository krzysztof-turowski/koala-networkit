#include <gtest/gtest.h>

#include <include/io/G6GraphReader.hpp>
#include <include/io/S6GraphReader.hpp>

#include <list>

struct GraphIOParameters {
    std::string G;
    int N;
    std::list<std::pair<int, int>> E;
};

class GraphReaderFromGraph6Test
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphReaderFromSparse6Test
    : public testing::TestWithParam<GraphIOParameters> { };

TEST_P(GraphReaderFromGraph6Test, test) {
  GraphIOParameters const& parameters = GetParam();
  Graph G = Koala::G6GraphReader(parameters.G);
  EXPECT_EQ(G.numberOfNodes(), parameters.N);
  EXPECT_EQ(G.numberOfEdges(), parameters.E.size());
  for (const auto &e : parameters.E) {
    EXPECT_TRUE(G.hasEdge(e.first, e.second));
  }
}

INSTANTIATE_TEST_CASE_P(
    test_small, GraphReaderFromGraph6Test, testing::Values(
        GraphIOParameters{
            "G?r@`_", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "G?qa`_", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 2}, {6, 1}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "GCQR@O", 8, {{3, 0}, {4, 1}, {5, 0}, {5, 3}, {6, 1}, {6, 2}, {7, 2}, {7, 4}}},
        GraphIOParameters{
            "Fw??G", 7, {{1, 0}, {2, 0}, {2, 1}, {6, 5}}}
));

TEST_P(GraphReaderFromSparse6Test, test) {
  GraphIOParameters const& parameters = GetParam();
  Graph G = Koala::S6GraphReader(parameters.G);
  EXPECT_EQ(G.numberOfNodes(), parameters.N);
  EXPECT_EQ(G.numberOfEdges(), parameters.E.size());
  for (const auto &e : parameters.E) {
    EXPECT_TRUE(G.hasEdge(e.first, e.second));
  }
}

INSTANTIATE_TEST_CASE_P(
    test_small, GraphReaderFromSparse6Test, testing::Values(
        GraphIOParameters{
            ":Fa@x^", 7, {{1, 0}, {2, 0}, {2, 1}, {6, 5}}},
        GraphIOParameters{
            ":Go@_YMb", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}}
));
