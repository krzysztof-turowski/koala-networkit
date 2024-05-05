#include <gtest/gtest.h>
#include <vector>

#include <max_cut/BranchAndBoundMaxCut.hpp>
#include <max_cut/RankTwoRelaxationMaxCut.hpp>
#include <max_cut/GoemansWilliamsonMaxCut.hpp>

struct MaxCutParameters {
    std::vector<std::vector<int>> graph;
    int expectedMaxCutValue;
};

class MaxCutTest : public testing::TestWithParam<MaxCutParameters> {};

TEST_P(MaxCutTest, TestMaxCutSolution) {
    MaxCutParameters const& parameters = GetParam();
    
    Koala::GoemansWilliamsonMaxCut solver(parameters.graph);
    
    solver.solve();
    
    EXPECT_EQ(solver.getMaxCutValue(), parameters.expectedMaxCutValue);
}

INSTANTIATE_TEST_SUITE_P(
    Default, MaxCutTest, testing::Values(
        // Test case 1
        MaxCutParameters{
            {
                {0, 2, 3, 0, 0},
                {2, 0, 4, 6, 0},
                {3, 4, 0, 5, 0},
                {0, 6, 5, 0, 1},
                {0, 0, 0, 1, 0}
            },
            17
        },
        // Test case 2
        MaxCutParameters{
            {
                {0, 1, 0, 1},
                {1, 0, 1, 0},
                {0, 1, 0, 1},
                {1, 0, 1, 0}
            },
            4
        },
        // Test case 3: Fully connected graph with uniform weights
        MaxCutParameters{
            {
                {0, 1, 1, 1},
                {1, 0, 1, 1},
                {1, 1, 0, 1},
                {1, 1, 1, 0}
            },
            4
        },

        // Test case 4: Star graph
        MaxCutParameters{
            {
                {0, 1, 1, 1, 1},
                {1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0}
            },
            4
        },

        // Test case 5: Linear chain graph
        MaxCutParameters{
            {
                {0, 1, 0, 0},
                {1, 0, 2, 0},
                {0, 2, 0, 3},
                {0, 0, 3, 0}
            },
            6
        },

        // Test case 6: Graph with varying edge weights
        MaxCutParameters{
            {
                {0, 10, 0, 0, 3},
                {10, 0, 2, 0, 0},
                {0, 2, 0, 4, 0},
                {0, 0, 4, 0, 5},
                {3, 0, 0, 5, 0}
            },
            22
        }
    )
);
