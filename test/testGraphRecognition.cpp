#include <gtest/gtest.h>

#include <list>

#include <recognition/PerfectGraphRecognition.hpp>

#include "helpers.hpp"

struct GraphRecognitionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    bool is_recognized;
};

class PerfectGraphRecognitionTest
    : public testing::TestWithParam<GraphRecognitionParameters> { };

TEST_P(PerfectGraphRecognitionTest, test) {
    GraphRecognitionParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::PerfectGraphRecognition(G);
    algorithm.run();

    auto is_perfect = algorithm.isPerfect();
    EXPECT_EQ(is_perfect, parameters.is_recognized);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, PerfectGraphRecognitionTest, testing::Values(
        GraphRecognitionParameters{4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, true},
        GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, false},
        GraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {3, 1}, {4, 3}}, true},
        GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {2, 4}, {4, 5}, {5, 6},{2,6},{2,3},{3,4},{0,2},{1,4}}, true}
));
