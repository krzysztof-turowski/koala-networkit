#include <gtest/gtest.h>

#include <list>

#include <flow/MaximumFlow.hpp>

#include "helpers.hpp"

struct MaximumFlowParameters {
    int N;
    std::list<std::tuple<int, int, int>> EW;
    int s, t;
    int flowSize;
};

class KingRaoTarjanMaximumFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(KingRaoTarjanMaximumFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW);
    auto algorithm = Koala::KingRaoTarjanMaximumFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, KingRaoTarjanMaximumFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23}
));
