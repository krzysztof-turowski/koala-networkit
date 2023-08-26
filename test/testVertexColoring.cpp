#include <gtest/gtest.h>

#include <list>

#include <coloring/ExactVertexColoring.hpp>
#include <coloring/GreedyVertexColoring.hpp>
#include <coloring/PerfectGraphVertexColoring.hpp>

#include "helpers.hpp"

struct VertexColoringParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int colors;
};

class RandomSequentialVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class LargestFirstVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class SmallestLastVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class SaturatedLargestFirstVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class GreedyIndependentSetVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class PerfectGraphVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class BrownEnumerationVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class ChristofidesEnumerationVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class BrelazEnumerationVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

class KormanEnumerationVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

auto test_set_exact = testing::Values(
    VertexColoringParameters{4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, 2},
    VertexColoringParameters{6, {{0, 1}, {0, 2}, {1, 2}, {0, 3}, {3, 4}, {1, 5}, {4, 5}}, 3},
    VertexColoringParameters{10,
        {{0, 2}, {0, 3}, {0, 4}, {0, 5}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 6}, {2, 8}, {3, 4},
            {3, 6}, {4, 6}, {5, 6}, {5, 7}, {5, 9}, {7, 8}, {7, 9}, {8, 9}},
        3},
    VertexColoringParameters{8,
        {{0, 1}, {0, 2}, {0, 4}, {0, 6}, {1, 2}, {1, 3}, {1, 7}, {2, 3}, {2, 4}, {3, 5}, {3, 7},
            {4, 5}, {4, 6}, {5, 6}, {5, 7}, {6, 7}},
        4},
    VertexColoringParameters{8,
        {{0, 1}, {0, 2}, {0, 4}, {0, 6}, {1, 3}, {1, 5}, {1, 7}, {2, 3}, {2, 4}, {2, 5}, {2, 6},
            {3, 4}, {3, 5}, {3, 7}, {4, 6}, {4, 7}, {5, 6}, {5, 7}, {6, 7}},
        4},
    VertexColoringParameters{10,
        {{0, 1}, {0, 2}, {0, 3}, {0, 5}, {0, 6}, {0, 7}, {1, 2}, {1, 3}, {1, 4}, {1, 6}, {1, 7},
            {2, 3}, {2, 4}, {2, 5}, {2, 7}, {3, 4}, {3, 5}, {3, 6}, {4, 5}, {4, 7}, {4, 8},
            {4, 9}, {5, 6}, {5, 8}, {5, 9}, {6, 7}, {6, 8}, {6, 9}, {7, 8}, {7, 9}, {8, 9}},
        5},
    VertexColoringParameters{4, {{0, 2}, {1, 3}, {2, 3}}, 2},
    VertexColoringParameters{9,
        {{0, 4}, {0, 5}, {0, 6}, {0, 8}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {2, 6}, {3, 7}, {3, 8},
            {4, 5}, {4, 7}, {4, 8}, {5, 7}, {5, 8}, {6, 7}, {7, 8}},
        4});

void check(const auto &parameters, const auto &colors) {
    for (auto [u, v] : parameters.E) {
        EXPECT_NE(colors.at(u), colors.at(v));
    }

    int max_color = 0;
    for (const auto &[v, c] : colors) {
        max_color = std::max(max_color, c);
    }
    EXPECT_EQ(max_color, parameters.colors);
}

TEST_P(RandomSequentialVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::RandomSequentialVertexColoring(G);
    algorithm.run();
    check(parameters, algorithm.getColoring());
}

INSTANTIATE_TEST_SUITE_P(
    test_example, RandomSequentialVertexColoringTest, testing::Values(
        VertexColoringParameters{4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, 2}
));

TEST_P(LargestFirstVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::LargestFirstVertexColoring(G);
    algorithm.run();
    check(parameters, algorithm.getColoring());
}

INSTANTIATE_TEST_SUITE_P(
    test_example, LargestFirstVertexColoringTest, testing::Values(
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

TEST_P(SmallestLastVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::SmallestLastVertexColoring(G);
    algorithm.run();
    check(parameters, algorithm.getColoring());
}

INSTANTIATE_TEST_SUITE_P(
    test_example, SmallestLastVertexColoringTest, testing::Values(
        VertexColoringParameters{
            8, {{0, 1}, {0, 2}, {0, 4}, {0, 6}, {1, 2}, {1, 3}, {1, 7}, {2, 3}, {2, 4}, {3, 5},
                {3, 7}, {4, 5}, {4, 6}, {5, 6}, {5, 7}, {6, 7}}, 5},
        VertexColoringParameters{
            8, {{0, 1}, {0, 2}, {0, 4}, {0, 6}, {1, 3}, {1, 5}, {1, 7}, {2, 3}, {2, 4}, {2, 5},
                {2, 6}, {3, 4}, {3, 5}, {3, 7}, {4, 6}, {4, 7}, {5, 6}, {5, 7}, {6, 7}}, 5},
        VertexColoringParameters{
            10, {{0, 2}, {0, 3}, {0, 4}, {0, 5}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 6}, {2, 8},
                 {3, 4}, {3, 6}, {4, 6}, {5, 6}, {5, 7}, {5, 9}, {7, 8}, {7, 9}, {8, 9}}, 4},
        VertexColoringParameters{
            10, {{0, 1}, {0, 2}, {0, 3}, {0, 5}, {0, 6}, {0, 7}, {1, 2}, {1, 3}, {1, 4}, {1, 6},
                 {1, 7}, {2, 3}, {2, 4}, {2, 5}, {2, 7}, {3, 4}, {3, 5}, {3, 6}, {4, 5}, {4, 7},
                 {4, 8}, {4, 9}, {5, 6}, {5, 8}, {5, 9}, {6, 7}, {6, 8}, {6, 9}, {7, 8}, {7, 9},
                 {8, 9}}, 6},
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

TEST_P(SaturatedLargestFirstVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::SaturatedLargestFirstVertexColoring(G);
    algorithm.run();
    check(parameters, algorithm.getColoring());
}

INSTANTIATE_TEST_SUITE_P(
    test_example, SaturatedLargestFirstVertexColoringTest, testing::Values(
        VertexColoringParameters{
            8, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 5}, {4, 6}, {4, 7}, {5, 6},
                {5, 7}, {6, 7}}, 4},
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

TEST_P(GreedyIndependentSetVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::GreedyIndependentSetVertexColoring(G);
    algorithm.run();
    check(parameters, algorithm.getColoring());
}

INSTANTIATE_TEST_SUITE_P(
    test_example, GreedyIndependentSetVertexColoringTest, testing::Values(
        VertexColoringParameters{
            6, {{0, 4}, {1, 4}, {2, 5}, {3, 5}, {4, 5}}, 3},
        VertexColoringParameters{
            10, {{0, 1}, {0, 4}, {0, 5}, {0, 6}, {0, 8}, {1, 2}, {1, 4}, {1, 5}, {1, 7}, {1, 8},
                 {2, 3}, {2, 4}, {2, 6}, {2, 7}, {2, 9}, {3, 5}, {3, 6}, {3, 7}, {3, 8}, {4, 8},
                 {5, 6}, {5, 7}, {5, 8}, {5, 9}, {6, 9}, {7, 9}, {8, 9}}, 5}
));

TEST_P(PerfectGraphVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::PerfectGraphVertexColoring(G);
    algorithm.run();
    check(parameters, algorithm.getColoring());
}

INSTANTIATE_TEST_SUITE_P(
    test_example, PerfectGraphVertexColoringTest, testing::Values(
        VertexColoringParameters{
            6, {{0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 4}, {2, 5}}, 3},
        VertexColoringParameters{
            6, {{0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 4}, {2, 5}, {3, 4}, {3, 5}, {4, 5}}, 3}
));

TEST_P(BrownEnumerationVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::BrownEnumerationVertexColoring(G);
    algorithm.run();

    auto colors = algorithm.getColoring();

    check(parameters, colors);
}

INSTANTIATE_TEST_SUITE_P(test_example,
    BrownEnumerationVertexColoringTest,
    test_set_exact);

TEST_P(ChristofidesEnumerationVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::ChristofidesEnumerationVertexColoring(G);
    algorithm.run();

    auto colors = algorithm.getColoring();

    check(parameters, colors);
}

INSTANTIATE_TEST_SUITE_P(test_example,
    ChristofidesEnumerationVertexColoringTest,
    test_set_exact);

TEST_P(BrelazEnumerationVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::BrelazEnumerationVertexColoring(G);
    algorithm.run();

    auto colors = algorithm.getColoring();

    check(parameters, colors);
}

INSTANTIATE_TEST_SUITE_P(test_example,
    BrelazEnumerationVertexColoringTest,
    test_set_exact);

TEST_P(KormanEnumerationVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::KormanEnumerationVertexColoring(G);
    algorithm.run();

    auto colors = algorithm.getColoring();

    check(parameters, colors);
}

INSTANTIATE_TEST_SUITE_P(test_example,
    KormanEnumerationVertexColoringTest,
    test_set_exact);
