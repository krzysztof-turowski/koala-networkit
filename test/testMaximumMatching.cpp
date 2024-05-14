#include <gtest/gtest.h>

#include <list>

#include <matching/MaximumMatching.hpp>

using Edge = std::pair<int, int>;
using WeightedEdge = std::tuple<int, int, int>;

enum MaxmimumWeightMatchingAlgorithm { edmonds, gabow, micali, scaling };

struct MaximumWeightMatchingParameters {
    NetworKit::count N;
    std::list<WeightedEdge> E;
    int maximumMatching;
};

struct MaximumCardinalityMatchingParameters {
    NetworKit::count N;
    std::list<Edge> E;
    NetworKit::count maximumMatching;
};

class MaximumWeightMatchingTest
    : public testing::TestWithParam<
        std::tuple<MaximumWeightMatchingParameters, MaxmimumWeightMatchingAlgorithm>> { };

class MicaliVaziraniTest
    : public testing::TestWithParam<MaximumCardinalityMatchingParameters> { };

NetworKit::Graph build_graph(const int &N, const std::list<Edge> &E) {
    NetworKit::Graph G(N, false, false);
    for (const auto &[u, v] : E) {
        G.addEdge(u, v);
    }
    return G;
}

NetworKit::Graph build_weighted_graph(const int &N, const std::list<WeightedEdge> &E) {
    NetworKit::Graph G(N, true, false);
    for (const auto &[u, v, w] : E) {
        G.addEdge(u, v, w);
    }
    return G;
}

void check_proper_matching(const auto& graph, const auto& matching) {
    for (auto [u, v] : matching) {
        if (v != NetworKit::none) {
            EXPECT_TRUE(graph.hasEdge(u, v));
            EXPECT_EQ(u, matching[v]);
        }
    }
}

int caluculate_matching_weight(const auto& graph, const auto& matching) {
    int weight = 0;
    for (auto [u, v] : matching) {
        if (v != NetworKit::none) {
            weight += static_cast<int>(graph.weight(u, v));
        }
    }
    return weight / 2;
}

NetworKit::count calculate_matching_size(const auto& matching) {
    NetworKit::count size = 0;
    for (auto [u, v] : matching) {
        if (v != NetworKit::none)
            size++;
    }
    return size / 2;
}

TEST_P(MaximumWeightMatchingTest, test) {
    auto const& [matchingParams, matchingAlgorithm] = GetParam();
    NetworKit::Graph G = build_weighted_graph(matchingParams.N, matchingParams.E);
    G.indexEdges(true);
    Koala::MaximumWeightMatching* algorithm = nullptr;

    switch (matchingAlgorithm) {
        case edmonds:
            algorithm = new Koala::EdmondsMaximumMatching(G);
            break;
        case gabow:
            algorithm = new Koala::GabowMaximumMatching(G);
            break;
        case micali:
            algorithm = new Koala::GalilMicaliGabowMaximumMatching(G);
            break;
        case scaling:
            algorithm = new Koala::GabowScalingMatching(G);
            break;
    }

    algorithm->run();
    auto matching = algorithm->getMatching();

    check_proper_matching(G, matching);

    auto matchingWeight = caluculate_matching_weight(G, matching);
    EXPECT_EQ(matchingParams.maximumMatching, matchingWeight);

    delete algorithm;
}

TEST_P(MicaliVaziraniTest, test) {
    auto const& matchingParams = GetParam();
    NetworKit::Graph G = build_graph(matchingParams.N, matchingParams.E);
    G.indexEdges(true);
    Koala::MicaliVaziraniMatching algorithm(G);

    algorithm.run();
    auto matching = algorithm.getMatching();

    check_proper_matching(G, matching);

    auto matching_size = calculate_matching_size(matching);
    EXPECT_EQ(matchingParams.maximumMatching, matching_size);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, MaximumWeightMatchingTest, testing::Combine(
        testing::Values(
            MaximumWeightMatchingParameters{3, {}, 0},
            MaximumWeightMatchingParameters{3, {{0, 1, 1}}, 1},
            MaximumWeightMatchingParameters{3, {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}}, 1},
            MaximumWeightMatchingParameters{
                12, {{1, 0, 4}, {2, 1, 2}, {4, 0, 4}, {4, 1, 12}, {4, 2, 3}, {6, 1, 17}, {6, 4, 6},
                     {6, 5, 19}, {7, 3, 13}, {8, 5, 5}, {8, 7, 12}, {9, 0, 16}, {9, 3, 11},
                     {10, 1, 8}, {10, 2, 3}, {10, 4, 14}, {10, 5, 3}, {10, 9, 17}, {11, 0, 12},
                     {11, 3, 16}, {11, 4, 3}, {11, 6, 19}, {11, 7, 1}, {11, 10, 13}}, 79},
            MaximumWeightMatchingParameters{
                5, {{1, 0, 6}, {3, 0, 5}, {3, 2, 3}, {4, 0, 8}, {4, 1, 5}}, 11},
            MaximumWeightMatchingParameters{
                4, {{2, 0, 8}, {3, 0, 6}, {3, 1, 10}, {3, 2, 10}}, 18},
            MaximumWeightMatchingParameters{
                12, {{3, 1, 37}, {4, 0, 73}, {4, 2, 53}, {4, 3, 34}, {5, 1, 24}, {5, 2, 16},
                     {6, 0, 1}, {6, 3, 72}, {6, 5, 4}, {7, 0, 91}, {7, 1, 47}, {7, 3, 44},
                     {7, 4, 72}, {7, 5, 51}, {8, 2, 87}, {8, 3, 21}, {8, 4, 92}, {8, 5, 54},
                     {8, 6, 67}, {9, 0, 31}, {9, 5, 71}, {9, 6, 49}, {10, 2, 85}, {10, 5, 63},
                     {10, 6, 90}, {11, 2, 38}, {11, 3, 64}, {11, 4, 40}, {11, 9, 11},
                     {11, 10, 73}}, 432},
            MaximumWeightMatchingParameters{
                7, {{1, 0, 7}, {2, 0, 3}, {2, 1, 4}, {3, 0, 9}, {3, 1, 8}, {4, 0, 6}, {5, 0, 2},
                     {5, 1, 7}, {5, 3, 4}, {6, 0, 10}, {6, 2, 8}, {6, 3, 10}, {6, 4, 4}}, 24},
            MaximumWeightMatchingParameters{
                6, {{1, 0, 6}, {2, 0, 7}, {2, 1, 10}, {3, 1, 2}, {4, 1, 9}, {5, 2, 8}}, 17},
            MaximumWeightMatchingParameters{
                10, {{2, 0, 12}, {3, 1, 3}, {5, 0, 28}, {6, 0, 1}, {6, 3, 15}, {7, 0, 3},
                     {7, 1, 12}, {8, 2, 13}, {8, 5, 30}, {8, 6, 16}, {8, 7, 15}, {9, 0, 23},
                     {9, 1, 15}, {9, 2, 28}, {9, 7, 20}}, 86},
            MaximumWeightMatchingParameters{
                5, {{1, 0, 7}, {2, 1, 8}, {3, 1, 5}, {3, 2, 9}, {4, 0, 3}, {4, 2, 4},
                    {4, 3, 7}}, 16},
            MaximumWeightMatchingParameters{
                9, {{1, 0, 11}, {3, 0, 4}, {3, 2, 20}, {4, 1, 15}, {4, 2, 11}, {5, 0, 8},
                    {5, 1, 16}, {5, 3, 20}, {5, 4, 15}, {6, 1, 19}, {6, 4, 13}, {6, 5, 18},
                    {7, 1, 18}, {7, 4, 18}, {8, 6, 10}}, 67},
            MaximumWeightMatchingParameters{
                6, {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {0, 4, 1}, {1, 2, 1}, {1, 3, 1}, {2, 3, 1},
                    {2, 5, 1}}, 3},
            MaximumWeightMatchingParameters{
                5, {{1, 0, 9}, {2, 1, 8}, {3, 0, 5}, {3, 1, 10}, {3, 2, 9}, {4, 0, 10},
                    {4, 3, 2}}, 20},
            MaximumWeightMatchingParameters{
                8, {{2, 0, 8}, {2, 1, 4}, {3, 0, 8}, {3, 1, 4}, {3, 2, 7}, {4, 1, 2}, {4, 2, 9},
                    {5, 2, 2}, {5, 3, 2}, {5, 4, 8}, {6, 0, 4}, {6, 1, 8}, {6, 2, 10}, {6, 3, 4},
                    {6, 4, 2}, {6, 5, 6}, {7, 1, 1}, {7, 2, 3}, {7, 6, 8}}, 28},
            MaximumWeightMatchingParameters{
                5, {{1, 0, 3}, {2, 0, 4}, {2, 1, 9}, {3, 0, 5}, {3, 2, 9}, {4, 1, 6}, {4, 2, 9},
                    {4, 3, 1}}, 15},
            MaximumWeightMatchingParameters{
                5, {{2, 0, 10}, {3, 0, 5}, {3, 1, 8}, {3, 2, 1}, {4, 0, 9}, {4, 1, 6}, {4, 2, 5},
                    {4, 3, 6}}, 18},
            MaximumWeightMatchingParameters{
                7, {{1, 0, 6}, {2, 1, 6}, {3, 2, 18}, {4, 0, 18}, {4, 1, 1}, {4, 2, 17}, {5, 0, 12},
                    {5, 2, 17}, {5, 3, 18}, {5, 4, 17}, {6, 0, 12}, {6, 3, 14}}, 49},
            MaximumWeightMatchingParameters{
                7, {{2, 0, 16}, {4, 0, 15}, {4, 2, 15}, {4, 3, 17}, {5, 1, 13}, {5, 2, 15},
                    {5, 3, 18}, {5, 4, 9}, {6, 0, 7}, {6, 1, 14}, {6, 3, 5}, {6, 5, 15}}, 48},
            MaximumWeightMatchingParameters{
                9, {{2, 1, 27}, {3, 0, 30}, {3, 2, 15}, {4, 0, 27}, {4, 3, 40}, {5, 3, 43},
                    {6, 0, 22}, {6, 1, 10}, {6, 5, 31}, {7, 0, 24}, {7, 4, 19}, {7, 6, 20},
                    {8, 2, 28}, {8, 4, 33}}, 127},
            MaximumWeightMatchingParameters{
                12, {{0, 1, 14}, {0, 4, 14}, {0, 9, 7}, {1, 3, 14}, {2, 3, 14}, {2, 4, 14},
                     {3, 5, 12}, {3, 6, 11}, {3, 8, 12}, {4, 7, 12}, {5, 6, 11}, {5, 11, 7},
                     {6, 8, 11}, {7, 8, 14}, {7, 10, 9}}, 62},
            MaximumWeightMatchingParameters{
                5, {{1, 0, 9}, {2, 0, 6}, {2, 1, 8}, {3, 0, 7}, {3, 1, 10}, {4, 0, 5}, {4, 2, 5},
                    {4, 3, 3}}, 16},
            MaximumWeightMatchingParameters{
                8, {{2, 1, 15}, {3, 1, 3}, {5, 2, 10}, {5, 4, 4}, {6, 0, 17}, {7, 0, 19}, {7, 1, 2},
                    {7, 2, 16}, {7, 6, 8}}, 40},
            MaximumWeightMatchingParameters{
                5, {{1, 0, 4}, {2, 0, 8}, {2, 1, 9}, {3, 0, 10}, {3, 1, 10}, {3, 2, 9}, {4, 0, 6},
                    {4, 1, 8}}, 19},
            MaximumWeightMatchingParameters{
                6, {{2, 0, 10}, {2, 1, 9}, {3, 2, 4}, {4, 0, 9}, {4, 3, 5}, {5, 1, 1}, {5, 3, 1},
                    {5, 4, 5}}, 19},
            MaximumWeightMatchingParameters{
                5, {{1, 0, 4}, {2, 0, 3}, {3, 1, 8}, {3, 2, 9}, {4, 0, 6}, {4, 1, 6}, {4, 2, 6},
                    {4, 3, 9}}, 15},
            MaximumWeightMatchingParameters{
                10, {{1, 0, 2}, {2, 0, 8}, {2, 1, 10}, {4, 1, 2}, {5, 0, 5}, {5, 1, 9}, {5, 2, 8},
                     {6, 0, 5}, {6, 2, 8}, {6, 4, 8}, {7, 2, 2}, {7, 3, 4}, {7, 4, 8}, {7, 5, 4},
                     {8, 5, 4}, {8, 7, 5}, {9, 0, 10}, {9, 1, 4}, {9, 7, 6}, {9, 8, 1}}, 36}),
        testing::Values(edmonds, gabow, micali, scaling)
));

INSTANTIATE_TEST_SUITE_P(
    test_example, MicaliVaziraniTest, testing::Values(
        MaximumCardinalityMatchingParameters{3, {}, 0},
        MaximumCardinalityMatchingParameters{3, {{0, 1}}, 1},
        MaximumCardinalityMatchingParameters{3, {{0, 1}, {1, 2}, {2, 0}}, 1},
        MaximumCardinalityMatchingParameters{
            12, {{3, 1}, {4, 0}, {4, 2}, {4, 3}, {5, 1}, {5, 2}, {6, 0}, {6, 3}, {6, 5}, {7, 0},
                 {7, 1}, {7, 3}, {7, 4}, {7, 5}, {8, 2}, {8, 3}, {8, 4}, {8, 5}, {8, 6}, {9, 0},
                 {9, 5}, {9, 6}, {10, 2}, {10, 5}, {10, 6}, {11, 2}, {11, 3}, {11, 4}, {11, 9},
                 {11, 10}}, 6},
        MaximumCardinalityMatchingParameters{
            14, {{9, 11}, {0, 12}, {8, 11}, {5, 10}, {10, 13}, {5, 7}, {1, 4}, {6, 7}, {3, 12},
                 {10, 12}, {0, 13}, {5, 13}, {2, 12}, {7, 11}, {6, 12}, {7, 8}, {0, 8}}, 6},
        MaximumCardinalityMatchingParameters{
            7, {{1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 0}, {5, 0}, {5, 1}, {5, 3}, {6, 0},
                {6, 2}, {6, 3}, {6, 4}}, 3},
        MaximumCardinalityMatchingParameters{
            14, {{1, 2}, {9, 11}, {0, 3}, {3, 4}, {6, 13}, {10, 12}, {0, 8}, {0, 13}, {4, 7},
                 {1, 10}, {0, 11}, {2, 11}, {10, 11}, {0, 5}, {6, 11}, {2, 8}, {7, 12}, {6, 8},
                 {0, 7}, {3, 7}}, 7},
        MaximumCardinalityMatchingParameters{
            14, {{4, 12}, {1, 11}, {6, 11}, {0, 10}, {0, 2}, {8, 13}, {0, 4}, {7, 12}, {1, 13},
                 {2, 4}, {6, 8}, {0, 5}, {1, 2}, {3, 7}}, 6},
        MaximumCardinalityMatchingParameters{
            6, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {2, 3}, {2, 5}}, 3},
        MaximumCardinalityMatchingParameters{
            8, {{2, 0}, {2, 1}, {3, 0}, {3, 1}, {3, 2}, {4, 1}, {4, 2}, {5, 2}, {5, 3}, {5, 4},
                {6, 0}, {6, 1}, {6, 2}, {6, 3}, {6, 4}, {6, 5}, {7, 1}, {7, 2}, {7, 6}}, 4},
        MaximumCardinalityMatchingParameters{
            12, {{0, 1}, {0, 4}, {0, 9}, {1, 3}, {2, 3}, {2, 4}, {3, 5}, {3, 6}, {3, 8}, {4, 7},
                 {5, 6}, {5, 11}, {6, 8}, {7, 8}, {7, 10}}, 6},
        MaximumCardinalityMatchingParameters{
            13, {{3, 1}, {4, 1}, {4, 2}, {5, 2}, {6, 2}, {6, 3}, {7, 2}, {7, 6}, {8, 1}, {8, 2},
                 {8, 4}, {8, 7}, {9, 1}, {9, 2}, {9, 3}, {9, 6}, {9, 7}, {10, 0}, {11, 1}, {11, 7},
                 {11, 8}, {12, 0}, {12, 5}, {12, 6}, {12, 7}}, 6},
        MaximumCardinalityMatchingParameters{
            15, {{8, 11}, {3, 5}, {7, 10}, {10, 13}, {6, 13}, {11, 13}, {12, 13}, {3, 4}, {1, 13},
                 {2, 11}, {0, 4}, {2, 3}, {2, 4}, {5, 9}, {0, 6}, {8, 13}, {7, 9}, {9, 14}, {3, 14},
                 {1, 2}, {9, 13}, {11, 14}, {1, 4}, {3, 13}, {4, 7}}, 7},
        MaximumCardinalityMatchingParameters{
            10, {{1, 0}, {2, 0}, {2, 1}, {4, 1}, {5, 0}, {5, 1}, {5, 2}, {6, 0}, {6, 2}, {6, 4},
                 {7, 2}, {7, 3}, {7, 4}, {7, 5}, {8, 5}, {8, 7}, {9, 0}, {9, 1}, {9, 7}, {9, 8}},
                 5},
        MaximumCardinalityMatchingParameters{
            13, {{3, 1}, {3, 2}, {5, 1}, {5, 3}, {7, 0}, {7, 4}, {8, 0}, {8, 5}, {8, 6}, {10, 0},
                 {10, 1}, {10, 2}, {10, 3}, {10, 7}, {10, 8}, {11, 0}, {11, 3}, {11, 8}, {12, 2},
                 {12, 10}}, 6},
        MaximumCardinalityMatchingParameters{
            8, {{2, 0}, {2, 1}, {3, 1}, {3, 2}, {4, 1}, {5, 0}, {5, 2}, {5, 3}, {6, 0}, {6, 1},
                {6, 2}, {6, 5}, {7, 1}, {7, 2}, {7, 5}, {7, 6}}, 4},
        MaximumCardinalityMatchingParameters{
            7, {{1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 1}, {4, 2}, {4, 3}, {5, 2}, {6, 0},
                {6, 1}, {6, 3}}, 3},
        MaximumCardinalityMatchingParameters{
            10, {{1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {3, 2}, {4, 1}, {4, 2}, {5, 0}, {5, 1},
                 {5, 2}, {5, 3}, {5, 4}, {6, 0}, {6, 1}, {6, 2}, {6, 3}, {6, 4}, {6, 5}, {7, 0},
                 {7, 1}, {7, 2}, {7, 3}, {7, 4}, {7, 5}, {7, 6}, {8, 0}, {8, 1}, {8, 3}, {8, 4},
                 {8, 5}, {8, 6}, {8, 7}, {9, 0}, {9, 1}, {9, 2}, {9, 4}, {9, 5}, {9, 6}, {9, 7},
                 {9, 8}}, 5},
        MaximumCardinalityMatchingParameters{
            14, {{1, 0}, {2, 0}, {3, 0}, {4, 0}, {4, 1}, {4, 2}, {5, 0}, {5, 2}, {5, 3}, {6, 3},
                 {7, 1}, {7, 2}, {7, 6}, {8, 2}, {8, 5}, {8, 6}, {9, 3}, {9, 4}, {10, 8}, {10, 9},
                 {11, 3}, {11, 8}, {12, 1}, {12, 5}, {12, 6}, {13, 1}, {13, 4}}, 7},
        MaximumCardinalityMatchingParameters{
            12, {{1, 0}, {2, 0}, {3, 0}, {3, 2}, {4, 0}, {4, 2}, {5, 1}, {5, 3}, {5, 4}, {6, 2},
                 {6, 3}, {7, 1}, {7, 2}, {7, 5}, {7, 6}, {8, 1}, {8, 2}, {8, 3}, {8, 5}, {9, 1},
                 {9, 2}, {9, 4}, {9, 5}, {9, 6}, {9, 8}, {10, 1}, {10, 3}, {10, 6}, {10, 7},
                 {10, 8}, {10, 9}, {11, 0}, {11, 1}, {11, 2}, {11, 3}, {11, 4}, {11, 5}, {11, 8},
                 {11, 9}}, 6},
        MaximumCardinalityMatchingParameters{
            13, {{1, 0}, {2, 0}, {3, 1}, {3, 2}, {4, 3}, {5, 3}, {6, 3}, {6, 4}, {7, 4}, {9, 1},
                 {9, 4}, {9, 5}, {11, 0}, {11, 6}, {11, 7}, {11, 8}, {12, 4}, {12, 5}, {12, 9}}, 6},
        MaximumCardinalityMatchingParameters{
            10, {{1, 0}, {2, 0}, {3, 0}, {3, 1}, {3, 2}, {4, 0}, {4, 1}, {4, 3}, {5, 4}, {6, 0},
                 {6, 1}, {6, 2}, {7, 0}, {7, 1}, {7, 2}, {7, 3}, {8, 6}, {8, 7}, {9, 2}, {9, 6},
                 {9, 7}, {9, 8}}, 5},
        MaximumCardinalityMatchingParameters{
            6, {{2, 0}, {2, 1}, {3, 0}, {3, 1}, {3, 2}, {4, 1}, {4, 2}, {5, 0}, {5, 1}, {5, 2},
                {5, 3}, {5, 4}}, 3},
        MaximumCardinalityMatchingParameters{
            12, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}, {1, 2}, {1, 3},
                 {1, 4}, {1, 5}, {1, 6}, {1, 7}, {2, 3}, {2, 4}, {2, 5}, {2, 6}, {2, 7}, {2, 9},
                 {3, 4}, {3, 5}, {3, 6}, {3, 7}, {4, 5}, {4, 6}, {4, 7}, {4, 10}, {5, 6}, {5, 7},
                 {6, 7}, {6, 11}}, 6}
));
