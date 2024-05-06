#include <gtest/gtest.h>
#include <vector>

#include <min_cut/PushRelabelMinCut.hpp>
#include <min_cut/HaoOrlinMinCut.hpp>
#include <min_cut/KargerMinCut.hpp>
#include <min_cut/KargerSteinMinCut.hpp>

struct MinCutParameters {
    std::vector<std::vector<int>> graph; // Adjacency matrix representation of the graph
    int expectedMinCutValue; // Expected minimum cut value
};

class MinCutTest : public testing::TestWithParam<MinCutParameters> {};

TEST_P(MinCutTest, TestMinCutSolution) {
    MinCutParameters const& parameters = GetParam();
    
    // Koala::HaoOrlinMinCut solver(parameters.graph.size(), 0,parameters.graph.size() - 1, parameters.graph);
    Koala::KargerMinCut solver(parameters.graph.size(), 3, parameters.graph);
    
    solver.solve();
    
    EXPECT_EQ(solver.getMinCut(), parameters.expectedMinCutValue);
}

INSTANTIATE_TEST_SUITE_P(
    Default, MinCutTest, testing::Values(
        // Test case 1
        MinCutParameters{
            {
                {0, 2, 3, 0, 0},
                {2, 0, 4, 6, 0},
                {3, 4, 0, 5, 0},
                {0, 6, 5, 0, 1},
                {0, 0, 0, 1, 0}
            },
            1
        },
        // Test case 2
        MinCutParameters{
            {
                {0, 1, 0, 1},
                {1, 0, 1, 0},
                {0, 1, 0, 1},
                {1, 0, 1, 0}
            },
            2
        },
        // Test case 3: Fully connected graph with uniform weights
        MinCutParameters{
            {
                {0, 1, 1, 1},
                {1, 0, 1, 1},
                {1, 1, 0, 1},
                {1, 1, 1, 0}
            },
            3
        },
        // Test case 4: Star graph
        MinCutParameters{
            {
                {0, 1, 1, 1, 1},
                {1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0}
            },
            1
        },
        // Test case 5: Linear chain graph
        MinCutParameters{
            {
                {0, 1, 0, 0},
                {1, 0, 2, 0},
                {0, 2, 0, 3},
                {0, 0, 3, 0}
            },
            1
        },
        // Test case 6: Graph with varying edge weights
        MinCutParameters{
            {
                {0, 10, 0, 0, 3},
                {10, 0, 2, 0, 0},
                {0, 2, 0, 4, 0},
                {0, 0, 4, 0, 5},
                {3, 0, 0, 5, 0}
            },
            5
        }
    )
);
