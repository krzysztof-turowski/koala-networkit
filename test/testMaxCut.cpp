#include <gtest/gtest.h>

#include <list>

#include <max_cut/BranchAndBoundMaxCut.hpp>
#include <max_cut/RankTwoRelaxationMaxCut.hpp>
#include <max_cut/GoemansWilliamsonMaxCut.hpp>
#include <max_cut/NaiveMaxCut.hpp>

#include "helpers.hpp"

struct MaxCutParameters {
    int N;
    std::list<std::tuple<int, int, int>> EW;
    int expectedMaxCutValue;
};

class MaxCutTest : public testing::TestWithParam<MaxCutParameters> {};

TEST_P(MaxCutTest, TestMaxCutSolution) {
    MaxCutParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, false);

    Koala::GoemansWilliamsonMaxCut algorithm(G);
    
    algorithm.run();
    
    EXPECT_EQ(algorithm.getMaxCutValue(), parameters.expectedMaxCutValue);
}

INSTANTIATE_TEST_SUITE_P(
    Default, MaxCutTest, testing::Values(
        // Test case 1
        MaxCutParameters{
            5,
            {
                {0, 1, 2},
                {0, 2, 3},
                {1, 2, 4},
                {1, 3, 6},
                {2, 3, 5},
                {3, 4, 1}
            },
            17
        },
        // Test case 2
        MaxCutParameters{
            4,
            {
                {0, 1, 1},
                {0, 3, 1},
                {1, 2, 1},
                {2, 3, 1}
            },
            4
        },
        // Test case 3: Fully connected graph with uniform weights
        MaxCutParameters{
            4,
            {
                {0, 1, 1},
                {0, 2, 1},
                {0, 3, 1},
                {1, 2, 1},
                {1, 3, 1},
                {2, 3, 1}
            },
            4
        },
        // Test case 4: Star graph
        MaxCutParameters{
            5,
            {
                {0, 1, 1},
                {0, 2, 1},
                {0, 3, 1},
                {0, 4, 1}
            },
            4
        },
        // Test case 5: Linear chain graph
        MaxCutParameters{
            4,
            {
                {0, 1, 1},
                {1, 2, 2},
                {2, 3, 3}
            },
            6
        },
        // Test case 6: Graph with varying edge weights
        MaxCutParameters{
            5,
            {
                {0, 1, 10},
                {0, 4, 3},
                {1, 2, 2},
                {2, 3, 4},
                {3, 4, 5}
            },
            22
        }
    )
);
