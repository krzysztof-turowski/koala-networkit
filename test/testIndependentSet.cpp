#include <gtest/gtest.h>

#include <list>

#include <independent_set/IndependentSet.hpp>

#include "helpers.hpp"

struct IndependentSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

template <typename Algorithm>
class SimpleGraphs : public testing::Test {
public:
    virtual void SetUp() {
    }

protected:
    void verify(IndependentSetParameters& parameters) {
        NetworKit::Graph G = build_graph(parameters.N, parameters.E);
        auto algorithm = Algorithm(G);
        algorithm.run();

        auto independentSet = algorithm.getIndependentSet();
        for (const auto &[u, v] : parameters.E) {
            EXPECT_FALSE(independentSet.contains(u) && independentSet.contains(v));
        }
        EXPECT_EQ(independentSet.size(), parameters.expectedSetSize);
    }
};

/*
the following graphs are from:
https://mathworld.wolfram.com/IndependentSet.html
*/
TYPED_TEST_CASE_P(SimpleGraphs);

TYPED_TEST_P(SimpleGraphs, WheelGraphW_8) {
    IndependentSetParameters parameters =
    {8, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 0},
        {7, 0}, {7, 1}, {7, 2}, {7, 3}, {7, 4}, {7, 5}, {7, 6},
        },
    3};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, UtilityGraphK_3_3) {
    IndependentSetParameters parameters =
    {6, {
        {0, 3}, {0, 4}, {0, 5},
        {1, 3}, {1, 4}, {1, 5},
        {2, 3}, {2, 4}, {2, 5}
        },
    3};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, PetersenGraph) {
    IndependentSetParameters parameters =
    {10, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
        {0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9},
        {5, 7}, {5, 8}, {6, 8}, {6, 9}, {7, 9}
        },
    4};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, FruchtGraph) {
    IndependentSetParameters parameters =
    {12, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 0},
        {0, 7}, {1, 8}, {2, 8}, {3, 9}, {4, 9}, {5, 10}, {6, 10},
        {7, 8}, {11, 7}, {11, 9}, {11, 10}
    },
    5};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, TwoK5) {
    IndependentSetParameters parameters =
    {10, {
        {0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4},
        {5,6}, {5,7}, {5,8}, {5,9}, {6,7}, {6,8}, {6,9}, {7,8}, {7,9}, {8,9},
    },
    2};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, TestingGraph) {
    IndependentSetParameters parameters =
    {6, {
        {1,2}, {0,3}, {1,3}, {2,3}, {0,4}, {2,4}, {0,5}, {1,5},
    },
    3};
    this->verify(parameters);
}
REGISTER_TYPED_TEST_CASE_P(SimpleGraphs, WheelGraphW_8, UtilityGraphK_3_3, PetersenGraph, FruchtGraph, TwoK5, TestingGraph);

using Algorithms = testing::Types<
    Koala::BruteForceIndependentSet,
    Koala::Mis1IndependentSet,
    Koala::Mis2IndependentSet,
    Koala::Mis3IndependentSet,
    Koala::Mis4IndependentSet,
    Koala::Mis5IndependentSet,
    Koala::MeasureAndConquerIndependentSet>;
INSTANTIATE_TYPED_TEST_CASE_P(IndependentSet, SimpleGraphs, Algorithms);
