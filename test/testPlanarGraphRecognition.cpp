#include <gtest/gtest.h>

#include <list>
#include <iostream>
#include <recognition/planar/HopcroftTarjan.hpp>

#include "helpers.hpp"

struct PlanarGraphRecognitionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    bool is_recognized;
};

class PlanarGraphRecognitionTest : public testing::TestWithParam<PlanarGraphRecognitionParameters> { };

TEST_P(PlanarGraphRecognitionTest, test) {
    PlanarGraphRecognitionParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::HopcroftTarjan(G, false);
    algorithm.run();
   
    auto is_planar = algorithm.isPlanar(); 
    std::cout << is_planar << "\n";
    EXPECT_EQ(is_planar, parameters.is_recognized);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, PlanarGraphRecognitionTest, testing::Values(
        PlanarGraphRecognitionParameters{4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, true},
        PlanarGraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}, true},
        PlanarGraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, false},
        PlanarGraphRecognitionParameters{6, {{0, 3}, {0, 4}, {0, 5}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}}, false}));
