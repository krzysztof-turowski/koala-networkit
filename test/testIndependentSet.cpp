#include <gtest/gtest.h>
#include <list>

#include <independentSet/SimpleIndependentSet.hpp>

struct IndependentSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

class BruteForceIndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };
class Mis1IndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };
class Mis2IndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };
class Mis3IndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };
class Mis4IndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };
class Mis5IndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };

NetworKit::Graph build_graph(const int &N, const std::list<std::pair<int, int>> &E) {
    NetworKit::Graph G(N, false, false);
    for (const auto &[u, v] : E) {
        G.addEdge(u, v);
    }
    return G;
}

namespace {
/*
simpleGraphs except {0, 1} are from:
https://mathworld.wolfram.com/IndependentSet.html
*/
std::vector<IndependentSetParameters> simpleGraphs {
    {1, {
        
        },
    1},
    {2, {
        {0, 1},
        },
    1},
    {3, {
        {0, 1}, {1, 2}
        },
    2},
    {6, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, 
        },
    3},
    // {26, {
    //     // 0 solo // = 1

    //     {1, 2}, {2, 3}, // = 2
    //     {4, 5}, {5, 6}, {6, 4}, // = 1

    //     {7, 8}, {8, 9}, {9, 10}, // = 2
    //     {11, 12}, {12, 13}, {13, 14}, {14, 11}, // = 2

    //     {15, 16}, {16, 17}, {17, 18}, {18, 19}, // = 3
    //     {20, 21}, {21, 22}, {22, 23}, {23, 24}, {24, 25}, // = 3

    //     //{, }, 
    //     },
    // 14},
    {8, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 0},
        {7, 0}, {7, 1}, {7, 2}, {7, 3}, {7, 4}, {7, 5}, {7, 6},
        },
    3},
    {6, {
        {0, 3}, {0, 4}, {0, 5}, 
        {1, 3}, {1, 4}, {1, 5},
        {2, 3}, {2, 4}, {2, 5}
        },
    3},
    {10, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
        {0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9},
        {5, 7}, {5, 8}, {6, 8}, {6, 9}, {7, 9}
        },
    4},
    {12, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 0},
        {0, 7}, {1, 8}, {2, 8}, {3, 9}, {4, 9}, {5, 10}, {6, 10},
        {7, 8}, {11, 7}, {11, 9}, {11, 10}   
    },
    5}
};
} /* unnamed namespace */

TEST_P(BruteForceIndependentSetTest, test) {
    IndependentSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);

    auto algorithm = Koala::BruteForceIndependentSet(G);
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

INSTANTIATE_TEST_SUITE_P(
    test_example, BruteForceIndependentSetTest, testing::ValuesIn(simpleGraphs)
);


TEST_P(Mis1IndependentSetTest, test) {
    IndependentSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);

    auto algorithm = Koala::Mis1IndependentSet(G);
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

INSTANTIATE_TEST_SUITE_P(
    test_example, Mis1IndependentSetTest, testing::ValuesIn(simpleGraphs)
);

TEST_P(Mis5IndependentSetTest, test) {
    IndependentSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);

    auto algorithm = Koala::Mis5IndependentSet(G);
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

INSTANTIATE_TEST_SUITE_P(
    test_example, Mis5IndependentSetTest, testing::ValuesIn(simpleGraphs)
);
