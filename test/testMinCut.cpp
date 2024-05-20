#include <gtest/gtest.h>

#include <list>

#include <min_cut/PushRelabelMinCut.hpp>
#include <min_cut/HaoOrlinMinCut.hpp>
#include <min_cut/KargerMinCut.hpp>
#include <min_cut/KargerSteinMinCut.hpp>
#include <min_cut/StoerWagnerMinCut.hpp>

#include <flow/MaximumFlow.hpp>

#include "helpers.hpp"

struct MinCutParameters {
    int N;
    std::list<std::tuple<int, int, int>> EW;
    int expectedMinCutValue;
};

class MinCutTest : public testing::TestWithParam<MinCutParameters> {};

TEST_P(MinCutTest, TestMinCutSolution) {
    MinCutParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, false);
    
    Koala::StoerWagnerMinCut algorithm(G);
    
    algorithm.run();
    
    EXPECT_EQ(algorithm.getMinCutValue(), parameters.expectedMinCutValue);
}

INSTANTIATE_TEST_SUITE_P(
    Default, MinCutTest, testing::Values(
        // Test case 1
        MinCutParameters{
            5,
            {
                {0, 1, 2},
                {0, 2, 3},
                {1, 2, 4},
                {1, 3, 6},
                {2, 3, 5},
                {3, 4, 1}
            },
            1
        },

        // Test case 2
        MinCutParameters{
            4,
            {
                {0, 1, 1},
                {0, 3, 1},
                {1, 2, 1},
                {2, 3, 1}
            },
            2
        },

        // Test case 3: Fully connected graph with uniform weights
        MinCutParameters{
            4,
            {
                {0, 1, 1},
                {0, 2, 1},
                {0, 3, 1},
                {1, 2, 1},
                {1, 3, 1},
                {2, 3, 1}
            },
            3
        },

        // Test case 4: Star graph
        MinCutParameters{
            5,
            {
                {0, 1, 1},
                {0, 2, 1},
                {0, 3, 1},
                {0, 4, 1}
            },
            1
        },

        // Test case 5: Linear chain graph
        MinCutParameters{
            4,
            {
                {0, 1, 1},
                {1, 2, 2},
                {2, 3, 3}
            },
            1
        },

        // Test case 6: Graph with varying edge weights
        MinCutParameters{
            5,
            {
                {0, 1, 10},
                {0, 4, 3},
                {1, 2, 2},
                {2, 3, 4},
                {3, 4, 5}
            },
            5
        }
    )
);
