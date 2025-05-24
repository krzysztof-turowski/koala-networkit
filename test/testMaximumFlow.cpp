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
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
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
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4,{{0, 1, 3},{0, 2, 5},{1, 2, 2},{2, 1, 3},{1, 3, 7},{2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {
                {0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}
            }, 0, 21, 0}
));

class PushRelabelMaximumFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(PushRelabelMaximumFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::PushRelabel(G, parameters.s, parameters.t); 
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, PushRelabelMaximumFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4,{{0, 1, 3},{0, 2, 5},{1, 2, 2},{2, 1, 3},{1, 3, 7},{2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {
                {0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}
            }, 0, 21, 0}
));

class MKMFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(MKMFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::MKMFlow(G, parameters.s, parameters.t); 
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, MKMFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4,{{0, 1, 3},{0, 2, 5},{1, 2, 2},{2, 1, 3},{1, 3, 7},{2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {
                {0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}
            }, 0, 21, 0}
));

class BKFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(BKFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::BKFlow(G, parameters.s, parameters.t); 
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, BKFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4,{{0, 1, 3},{0, 2, 5},{1, 2, 2},{2, 1, 3},{1, 3, 7},{2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {
                {0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}
            }, 0, 21, 0}
));