#include <gtest/gtest.h>

#include <iostream>
#include <list>
#include <recognition/planar/PlanarGraphRecognition.hpp>

#include "helpers.hpp"

struct PlanarGraphRecognitionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    Koala::PlanarGraphRecognition::State is_recognized;
};

class PlanarGraphRecognitionTest
        : public testing::TestWithParam<PlanarGraphRecognitionParameters> {};

TEST_P(PlanarGraphRecognitionTest, test) {
    PlanarGraphRecognitionParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::HopcroftTarjan(G, false);
    algorithm.run();

    auto is_planar = algorithm.isPlanar();
    EXPECT_EQ(is_planar, parameters.is_recognized);
}

INSTANTIATE_TEST_SUITE_P(test_example, PlanarGraphRecognitionTest,
        testing::Values(
        PlanarGraphRecognitionParameters{
                4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, Koala::PlanarGraphRecognition::State::PLANAR},
        PlanarGraphRecognitionParameters{
                5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}, Koala::PlanarGraphRecognition::State::PLANAR},
        PlanarGraphRecognitionParameters{5,
                                         {{0, 1},
                                          {0, 2},
                                          {0, 3},
                                          {0, 4},
                                          {1, 2},
                                          {1, 3},
                                          {1, 4},
                                          {2, 3},
                                          {2, 4},
                                          {3, 4}},
                                         Koala::PlanarGraphRecognition::State::NOT_PLANAR}, // K5
        PlanarGraphRecognitionParameters{6,
                                         {{0, 3},
                                          {0, 4},
                                          {0, 5},
                                          {1, 3},
                                          {1, 4},
                                          {1, 5},
                                          {2, 3},
                                          {2, 4},
                                          {2, 5}},
                                         Koala::PlanarGraphRecognition::State::NOT_PLANAR})); // K3,3
