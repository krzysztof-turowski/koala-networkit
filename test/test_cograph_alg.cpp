//
// Created by scales on 19.03.24.
//

#include <gtest/gtest.h>
#include <list>

#include "cograph_rec/CographAlg.hpp"
#include "helpers.hpp"

struct GraphRecognitionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    bool is_cograph;
};

class cograph_alg_test
        : public testing::TestWithParam<GraphRecognitionParameters> {
};

TEST_P(cograph_alg_test, test) {
    GraphRecognitionParameters const &parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::CographRecognition(G);
    algorithm.run();

    auto is_сograph = algorithm.IsCograph();
    EXPECT_EQ(is_сograph, parameters.is_cograph);
}

INSTANTIATE_TEST_SUITE_P
(
        test_example, cograph_alg_test, testing::Values
        (
                GraphRecognitionParameters{3, {{0, 1}, {0, 2}}, true},
                GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {2, 4}}, false},
                GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {3, 4}}, true},
                GraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {1, 3}, {3, 4}}, false},
                GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {1, 4}, {2, 4}}, false},
                GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {0, 2}, {1, 4}, {2, 4}, {3, 2}, {3, 4},
                                               {4, 5}, {5, 6}, {2, 6}}, false},
                GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {0, 2}, {1, 4}, {2, 4}, {3, 2}, {3, 4},
                                               {4, 5}, {5, 6}}, false}
        )
);