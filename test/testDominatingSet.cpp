#include <gtest/gtest.h>

#include <dominating_set/ExactDominatingSet.hpp>
#include <set_cover/BranchAndReduceSetCover.hpp>

#include "helpers.hpp"

struct DominatingSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int minimumDominatingSetSize;
};

class GrandoniTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class FominGrandoniKratschTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class RooijBodlaenderTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class FominKratschWoegingerTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class SchiermeyerTest
    : public testing::TestWithParam<DominatingSetParameters> {};

auto parameter_set = DominatingSetParameters{
    32,
    {{0, 1}, {0, 2}, {0, 4}, {0, 8}, {0, 16}, {1, 2}, {1, 3}, {1, 5}, {1, 9}, {1, 17},
        {2, 3}, {2, 4}, {2, 6}, {2, 10}, {2, 18}, {3, 4}, {3, 5}, {3, 7}, {3, 11}, {3, 19},
        {4, 5}, {4, 6}, {4, 8}, {4, 12}, {4, 20}, {5, 6}, {5, 7}, {5, 9}, {5, 13}, {5, 21},
        {6, 7}, {6, 8}, {6, 10}, {6, 14}, {6, 22}, {7, 8}, {7, 9}, {7, 11}, {7, 15}, {7, 23},
        {8, 9}, {8, 10}, {8, 12}, {8, 16}, {8, 24}, {9, 10}, {9, 11}, {9, 13}, {9, 17}, {9, 25},
        {10, 11}, {10, 12}, {10, 14}, {10, 18}, {10, 26}, {11, 12}, {11, 13}, {11, 15},
        {11, 19}, {11, 27}, {12, 13}, {12, 14}, {12, 16}, {12, 20}, {12, 28}, {13, 14},
        {13, 15}, {13, 17}, {13, 21}, {13, 29}, {14, 15}, {14, 16}, {14, 18}, {14, 22},
        {14, 30}, {15, 16}, {15, 17}, {15, 19}, {15, 23}, {15, 31}, {16, 17}, {16, 18},
        {16, 20}, {16, 24}, {17, 18}, {17, 19}, {17, 21}, {17, 25}, {18, 19}, {18, 20},
        {18, 22}, {18, 26}, {19, 20}, {19, 21}, {19, 23}, {19, 27}, {20, 21}, {20, 22},
        {20, 24}, {20, 28}, {21, 22}, {21, 23}, {21, 25}, {21, 29}, {22, 23}, {22, 24},
        {22, 26}, {22, 30}, {23, 24}, {23, 25}, {23, 27}, {23, 31}, {24, 25}, {24, 26},
        {24, 28}, {25, 26}, {25, 27}, {25, 29}, {26, 27}, {26, 28}, {26, 30}, {27, 28},
        {27, 29}, {27, 31}, {28, 29}, {28, 30}, {29, 30}, {29, 31}, {30, 31}},
    5};

TEST_P(GrandoniTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, GrandoniTest, testing::Values(parameter_set));

TEST_P(FominGrandoniKratschTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, FominGrandoniKratschTest, testing::Values(parameter_set));

TEST_P(RooijBodlaenderTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, RooijBodlaenderTest, testing::Values(parameter_set));

TEST_P(FominKratschWoegingerTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::FominKratschWoegingerDominatingSet(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, FominKratschWoegingerTest, testing::Values(parameter_set));

TEST_P(SchiermeyerTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::SchiermeyerDominatingSet(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, SchiermeyerTest, testing::Values(parameter_set));
