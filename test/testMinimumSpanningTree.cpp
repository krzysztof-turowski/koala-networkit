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

    void test_approximate_mst() {
        SpanningTreeParameters const& parameters = GetParam();
        auto G = build_graph(parameters.N, parameters.EW, false);
        auto mst = Algorithm(G);
        
        int max_w = 0;
        G.forEdges([&max_w](NetworKit::node, NetworKit::node, NetworKit::edgeweight ew, NetworKit::edgeid) {
            max_w = std::max(max_w, static_cast<int>(ew));
        });
        
        for (float eps : std::array{0.005f, 0.01f, 0.025f, 0.05f, 0.15f, 0.25f, 0.35f, 0.45f}){
            mst.run(max_w, eps);
            float tree_weight = mst.getTreeWeight();

            // Logging information as this is an approxmiate algorithm
            // and the tests may be not pass in rare cases.
            std::cout << "Approx MST with eps = " << eps << std::endl;
            std::cout << "expected tree weight: " << parameters.treeWeight << std::endl;
            std::cout << "algorithm got: " << tree_weight << std::endl;
            std::cout << "-------------------------------" << std::endl;
            
            EXPECT_TRUE(parameters.treeWeight * (1 - eps) <= tree_weight);
            EXPECT_TRUE(parameters.treeWeight * (1 + eps) >= tree_weight);
        }
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

class ChazelleRubinfeldTrevisanMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree> { };

TEST_P(ChazelleRubinfeldTrevisanMinimumSpanningTreeTest, test_example) {
    test_approximate_mst();
}

INSTANTIATE_TEST_SUITE_P(test_example, ChazelleRubinfeldTrevisanMinimumSpanningTreeTest, example_trees);

class Chazelle2000MinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::Chazelle2000MinimumSpanningTree> { };

TEST_P(Chazelle2000MinimumSpanningTreeTest, test_example) {
    test_mst();
}

INSTANTIATE_TEST_CASE_P(test_example, Chazelle2000MinimumSpanningTreeTest, example_trees);