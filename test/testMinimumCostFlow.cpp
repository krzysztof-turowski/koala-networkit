#include <gtest/gtest.h>

#include <list>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <flow/minimum_cost_flow/SuccessiveApproximationMinimumCostFlow.hpp>
#include <flow/minimum_cost_flow/EdmondsKarpMinimumCostFlow.hpp>
#include <flow/minimum_cost_flow/OrlinMinimumCostFlow.hpp>

#include "test/helpers.hpp"

struct MinCostFlowParams {
    int N;
    // from, to, capacity, cost
    std::list<std::tuple<int, int, int, int>> EW;
    std::unordered_map<int, int> excess;
    int minCost;
};

Koala::MCFlowNetwork getInstance(MinCostFlowParams const& params) {
    std::list<std::tuple<int, int, int>> edges;

    std::unordered_map<NetworKit::Edge, std::int64_t> costs;
    std::unordered_map<NetworKit::node, std::int64_t> excess(
        params.excess.begin(), params.excess.end());

    for (auto [u, v, capacity, cost] : params.EW) {
        costs[{static_cast<NetworKit::node>(u), static_cast<NetworKit::node>(v)}] = cost;
        edges.push_back({u, v, capacity});
    }
    NetworKit::Graph G = build_graph(params.N, edges, true, false);

    return Koala::MCFlowNetwork(G, costs, excess);
}

const std::vector<MinCostFlowParams> basic_tests = {
    MinCostFlowParams{
        4,
        {{0, 2, 2, 1}, {2, 0, 1, 0}, {0, 3, 3, 1}, {2, 3, 2, 0}, {1, 2, 2, 1}, {1, 3, 2, 1}},
        {{0, 3}, {1, 2}, {3, -5}}, 5
    },
    MinCostFlowParams{
        8,
        {{1, 0, 3, 5}, {2, 0, 2, 1}, {0, 3, 6, 1}, {5, 4, 0, 0},
         {6, 4, 4, 3}, {6, 7, 1, 2}, {7, 4, 2, 1}},
        {{7, 2}, {6, 2}, {4, -4}, {0, -1}, {1, 2}, {2, 2}, {3, -3}}, 23
    },
    MinCostFlowParams{
        3,
        {{0, 1, 5, 1}, {1, 2, 5, 2}},
        {{0, 3}, {2, -3}},
        9
    },
    MinCostFlowParams{
        4,
        {{0, 1, 2, 1}, {1, 3, 2, 1}, {0, 2, 5, 5}, {2, 3, 5, 5}},
        {{0, 4}, {3, -4}},
        24
    },
    MinCostFlowParams{
        2,
        {{0, 1, 10, 3}},
        {{0, 2}, {1, -2}},
        6
    },
    MinCostFlowParams{
        5,
        {{0, 2, 3, 2}, {1, 2, 4, 1}, {2, 3, 5, 3}, {2, 4, 5, 4}},
        {{0, 2}, {1, 3}, {3, -1}, {4, -4}},
        26
    },
    MinCostFlowParams{
        4,
        {{0, 1, 5, 0}, {1, 2, 5, 0}, {2, 3, 5, 0}, {0, 3, 1, 10}},
        {{0, 4}, {3, -4}},
        0
    },
    MinCostFlowParams{
        4,
        {{0, 1, 10, 1}, {1, 2, 5, 5}, {2, 3, 10, 1}},
        {{0, 5}, {3, -5}},
        35
    }
};

const std::vector<MinCostFlowParams> complex_tests = {
    MinCostFlowParams{
        6,
        {
            {0, 3, 1, 10}, {0, 4, 1, 2}, {0, 5, 1, 8},
            {1, 3, 1, 9},  {1, 4, 1, 8}, {1, 5, 1, 1},
            {2, 3, 1, 2},  {2, 4, 1, 9}, {2, 5, 1, 8}
        },
        {{0, 1}, {1, 1}, {2, 1}, {3, -1}, {4, -1}, {5, -1}},
        5
    },
    MinCostFlowParams{
        4,
        {
            {0, 1, 1, 0}, {0, 2, 1, 2},
            {1, 2, 1, 1}, {1, 3, 1, 4},
            {2, 3, 1, 0}
        },
        {{0, 2}, {3, -2}},
        6
    },
    MinCostFlowParams{
        6,
        {
            {0, 1, 10, 2}, {0, 2, 5, 5},
            {1, 2, 5, 1},  {1, 3, 5, 4},
            {2, 3, 3, 1},  {2, 4, 8, 2},
            {3, 4, 5, 1},  {3, 5, 10, 3},
            {4, 5, 10, 2}
        },
        {{0, 12}, {5, -12}},
        98
    },
    MinCostFlowParams{
        5,
        {
            {0, 1, 10, 100}, {0, 2, 10, 1},
            {2, 3, 10, 1},   {3, 1, 10, 1},
            {1, 4, 10, 1}
        },
        {{0, 5}, {4, -5}},
        20
    },
    MinCostFlowParams{
        7,
        {
            {0, 3, 15, 2}, {1, 3, 10, 3}, {2, 3, 10, 1},
            {0, 4, 5, 10},
            {3, 4, 10, 2}, {3, 5, 10, 4}, {3, 6, 10, 3},
            {1, 5, 5, 6},  {2, 6, 5, 2}
        },
        {{0, 10}, {1, 5}, {2, 5}, {4, -8}, {5, -7}, {6, -5}},
        84
    }
};

class EdmondsKarpTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(EdmondsKarpTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto network = getInstance(parameters);
    auto algorithm = Koala::EdmondsKarpMinimumCostFlow(network);
    algorithm.run();

    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example_edmonds, EdmondsKarpTest, testing::ValuesIn(basic_tests));
INSTANTIATE_TEST_SUITE_P(test_complex_edmonds, EdmondsKarpTest, testing::ValuesIn(complex_tests));

class OrlinTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(OrlinTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto network = getInstance(parameters);
    auto algorithm = Koala::OrlinMinimumCostFlow(network);
    algorithm.run();

    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example_orlin, OrlinTest, testing::ValuesIn(basic_tests));
INSTANTIATE_TEST_SUITE_P(test_complex_orlin, OrlinTest, testing::ValuesIn(complex_tests));

class SuccessiveApproximationTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(SuccessiveApproximationTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto network = getInstance(parameters);
    auto algorithm = Koala::SuccessiveApproximationMinimumCostFlow(network);
    algorithm.run();

    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(
    test_example_successive, SuccessiveApproximationTest, testing::ValuesIn(basic_tests));
INSTANTIATE_TEST_SUITE_P(
    test_complex_successive, SuccessiveApproximationTest, testing::ValuesIn(complex_tests));
