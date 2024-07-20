#include <gtest/gtest.h>

#include <list>

#include <recognition/CographRecognitionOther.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

#include "helpers.hpp"

struct GraphRecognitionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    bool is_recognized;
};

auto perfect_graphs = testing::Values(
    GraphRecognitionParameters{4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, false},
    GraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {3, 1}, {4, 3}}, true},
    GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {2, 4}, {4, 5}, {5, 6}, {2, 6}, {2, 3}, {3, 4},
                                   {0, 2}, {1, 4}}, true}
);

class PerfectGraphRecognitionTest
    : public testing::TestWithParam<GraphRecognitionParameters> { };

TEST_P(PerfectGraphRecognitionTest, test_perfect_graphs) {
    auto parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::PerfectGraphRecognition(G);
    algorithm.run();

    auto is_perfect = algorithm.isPerfect();
    EXPECT_EQ(is_perfect, parameters.is_recognized);
}

INSTANTIATE_TEST_SUITE_P(
    test_perfect_graphs_example, PerfectGraphRecognitionTest, perfect_graphs);

auto cographs = testing::Values(
    GraphRecognitionParameters{4, {{0, 2}, {0, 3}, {1, 3}, {2, 3}}, true},
    GraphRecognitionParameters{5, {{0, 3}, {0, 4}, {3, 4}, {1, 3}, {2, 4}}, false},
    GraphRecognitionParameters{3, {{0, 1}, {0, 2}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {2, 4}}, false},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {3, 4}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {1, 3}, {3, 4}}, false},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {1, 4}, {2, 4}}, false},
    GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {0, 2}, {1, 4}, {2, 4}, {3, 2}, {3, 4},
                                   {4, 5}, {5, 6}, {2, 6}}, false},
    GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {0, 2}, {1, 4}, {2, 4}, {3, 2}, {3, 4},
                                   {4, 5}, {5, 6}}, false}
);

template <typename Algorithm>
void CographRecognitionTest(GraphRecognitionParameters parameters) {
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Algorithm(G);
    algorithm.run();

    auto is_cograph = algorithm.isCograph();
    EXPECT_EQ(is_cograph, parameters.is_recognized);
}

class CographHabibPaulRecognitionTest
    : public testing::TestWithParam<GraphRecognitionParameters> { };

TEST_P(CographBCHPRecognitionTest, test_cographs) {
    CographRecognitionTest<Koala::HabibPaulCographRecognition>(GetParam());
}

INSTANTIATE_TEST_SUITE_P(
    test_cographs_habib_paul_example, CographHabibPaulRecognitionTest, cographs);
