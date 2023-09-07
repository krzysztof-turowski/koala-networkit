#include <gtest/gtest.h>

#include <list>
#include <iostream>

#include <networkit/graph/GraphTools.hpp>

#include <mst/MinimumSpanningTree.hpp>

#include "helpers.hpp"

struct SpanningTreeParameters {
    int N;
    std::list<std::tuple<int, int, int>> EW;
    int treeWeight;
};

template <class Algorithm>
class MinimumSpanningTreeTest : public testing::TestWithParam<SpanningTreeParameters> {
 public:
    void test_mst() {
        SpanningTreeParameters const& parameters = GetParam();
        auto G = build_graph(parameters.N, parameters.EW, false);
        auto mst = Algorithm(G);
        mst.run();
        int tree_weight = 0;
        mst.getForest().forEdges(
            [&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
                tree_weight += G.weight(u, v);
            });
        EXPECT_EQ(parameters.N - 1, mst.getForest().numberOfEdges());
        EXPECT_EQ(tree_weight, parameters.treeWeight);
    }
};

auto example_trees = testing::Values(
    SpanningTreeParameters{
        4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 20}
);

class KruskalMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::KruskalMinimumSpanningTree> { };

TEST_P(KruskalMinimumSpanningTreeTest, test_example) {
    test_mst();
}

INSTANTIATE_TEST_SUITE_P(test_example, KruskalMinimumSpanningTreeTest, example_trees);

class BoruvkaMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::BoruvkaMinimumSpanningTree> { };

TEST_P(BoruvkaMinimumSpanningTreeTest, test_example) {
    test_mst();
}

INSTANTIATE_TEST_SUITE_P(test_example, BoruvkaMinimumSpanningTreeTest, example_trees);

class KargerKleinTarjanMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::KargerKleinTarjanMinimumSpanningTree> { };

TEST_P(KargerKleinTarjanMinimumSpanningTreeTest, test_example) {
    test_mst();
}

INSTANTIATE_TEST_SUITE_P(test_example, KargerKleinTarjanMinimumSpanningTreeTest, example_trees);
