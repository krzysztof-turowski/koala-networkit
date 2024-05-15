//
// Created by milana on 17.11.23.
//
#include <gtest/gtest.h>


#include <list>

#include "recognition/CographRecognition.hpp"

#include "helpers.hpp"

struct GraphRecognitionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    bool is_recognized;
};

class CographRecognitionTest
        : public testing::TestWithParam<GraphRecognitionParameters> {
};

class CographLexBfsRecognitionTest
        : public testing::TestWithParam<GraphRecognitionParameters> {
};

TEST_P(CographRecognitionTest, test) {
    GraphRecognitionParameters const &parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::CorneilStewartPerlCographRecognition(G);
    algorithm.run();

    auto is_cograph = algorithm.isCograph();
    EXPECT_EQ(is_cograph, parameters.is_recognized);
}

INSTANTIATE_TEST_SUITE_P(
        test_example, CographRecognitionTest, testing::Values(
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
));

TEST_P(CographLexBfsRecognitionTest, test) {
    GraphRecognitionParameters const &parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::BretscherCorneilHabibPaulCographRecognition(G);
    algorithm.run();

    auto is_cograph = algorithm.isCograph();
    EXPECT_EQ(is_cograph, parameters.is_recognized);
}

INSTANTIATE_TEST_SUITE_P(
        test_example, CographLexBfsRecognitionTest, testing::Values(
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
));
