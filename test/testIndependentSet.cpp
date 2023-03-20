#include <gtest/gtest.h>
#include <list>

#include <independentSet/SimpleIndependentSet.hpp>

struct IndependentSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

template <typename T> 
class SimpleGraphs : public testing::Test {
public:
    virtual void SetUp() {
    }

protected:
    void verify(IndependentSetParameters& parameters) {
        NetworKit::Graph G = build_graph(parameters.N, parameters.E);
        auto algorithm = T(G);
        algorithm.run();

        auto independentSet = algorithm.getIndependentSet();
        for (const auto &[u, v] : parameters.E) {
            EXPECT_FALSE(independentSet[u] && independentSet[v]);
        }

        int setSize = 0;
        for (const auto &[v, belongs] : independentSet) {
            setSize += (belongs);
        }
        EXPECT_EQ(setSize, parameters.expectedSetSize);
    }

    NetworKit::Graph build_graph(const int &N, const std::list<std::pair<int, int>> &E) {
        NetworKit::Graph G(N, false, false);
        for (const auto &[u, v] : E) {
            G.addEdge(u, v);
        }
        return G;
    }
};

TYPED_TEST_CASE_P(SimpleGraphs);

/*
the following graphs are from:
https://mathworld.wolfram.com/IndependentSet.html
*/

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

REGISTER_TYPED_TEST_CASE_P(SimpleGraphs, WheelGraphW_8, UtilityGraphK_3_3, PetersenGraph, FruchtGraph);

typedef testing::Types<
    Koala::BruteForceIndependentSet, 
    Koala::Mis1IndependentSet, 
    Koala::Mis2IndependentSet,
    Koala::Mis3IndependentSet,
    Koala::Mis4IndependentSet,
    Koala::Mis5IndependentSet> Algorithms;
INSTANTIATE_TYPED_TEST_CASE_P(IndependentSet, SimpleGraphs, Algorithms);
