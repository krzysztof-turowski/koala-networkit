#include <gtest/gtest.h>

#include <algorithm>
#include <list>
#include <iostream>
#include <stdexcept>

#include <networkit/graph/GraphTools.hpp>
#include "graph/GraphTools.hpp"
#include "mst/MinimumSpanningTree.hpp"

#include "test/helpers.hpp"

struct SpanningTreeParameters {
    int N;
    std::list<std::tuple<int, int, int>> EW;
    int treeWeight;
};

auto example_graphs = testing::Values(
    SpanningTreeParameters { 1, {}, 0 },
    SpanningTreeParameters { 2, { { 0, 1, 7 } }, 7 },
    SpanningTreeParameters { 3, { { 0, 1, 9 }, { 0, 2, 3 }, { 1, 2, 4 } }, 7 },
    SpanningTreeParameters {
        4, { { 0, 1, 1 }, { 1, 2, 1 }, { 2, 3, 1 }, { 3, 0, 1 }, { 0, 2, 5 } }, 3 },
    SpanningTreeParameters {
        4, { { 0, 1, 10 }, { 0, 2, 5 }, { 1, 2, 15 }, { 1, 3, 5 }, { 2, 3, 10 } }, 20 },
    SpanningTreeParameters {
        5, { { 0, 1, 8 }, { 0, 2, 2 }, { 1, 2, 3 }, { 1, 3, 1 }, { 2, 4, 4 },
             { 3, 4, 7 } }, 10 },
    SpanningTreeParameters {
        10, { {0, 3, 1}, {7, 8, 2}, {6, 9, 3}, {4, 5, 4}, {1, 2, 5}, {0, 1, 6},
              {2, 6, 7}, {8, 9, 8}, {5, 6, 9}, {4, 8, 10}, {1, 4, 11}, {3, 7, 12},
              {3, 4, 13} }, 45 },
    SpanningTreeParameters {
        17, { {0, 1, 81}, {0, 5, 92}, {0, 14, 167}, {0, 16, 1}, {1, 15, 3}, {1, 6, 14},
              {1, 11, 52}, {2, 6, 115}, {2, 14, 4}, {3, 1, 140}, {3, 4, 120}, {3, 6, 94},
              {3, 13, 5}, {4, 6, 32}, {4, 12, 6}, {5, 1, 181}, {5, 4, 165}, {5, 11, 7},
              {6, 10, 2}, {7, 0, 90}, {7, 2, 84}, {7, 9, 8}, {7, 11, 26}, {8, 5, 127},
              {8, 12, 193}, {8, 14, 9}, {9, 6, 102}, {9, 11, 113}, {9, 13, 103}, {10, 11, 122},
              {10, 14, 153}, {11, 6, 179}, {11, 15, 136}, {12, 13, 24}, {13, 1, 190},
              {13, 4, 148}, {13, 6, 44}, {14, 4, 75}, {14, 7, 64}, {14, 11, 93}, {15, 6, 192},
              {15, 12, 34}, {16, 2, 95} }, 338 },
    SpanningTreeParameters{
        17, { {0, 2, 509}, {0, 4, 574}, {0, 6, 478}, {0, 8, 31}, {0, 11, 293}, {0, 12, 977},
              {0, 16, 1}, {1, 0, 404}, {1, 3, 306}, {1, 5, 292}, {1, 9, 872}, {1, 10, 222},
              {1, 11, 980}, {1, 15, 2}, {2, 3, 656}, {2, 6, 792}, {2, 9, 743}, {2, 14, 3},
              {2, 15, 333}, {3, 4, 561}, {3, 5, 286}, {3, 6, 318}, {3, 7, 469}, {3, 8, 544},
              {3, 12, 847}, {3, 13, 4}, {4, 5, 798}, {4, 7, 936}, {4, 10, 933}, {4, 12, 5},
              {4, 14, 226}, {5, 6, 747}, {5, 8, 429}, {5, 11, 6}, {5, 12, 620}, {6, 7, 322},
              {6, 10, 7}, {7, 8, 756}, {7, 9, 8}, {7, 11, 776}, {8, 14, 30}, {9, 12, 617},
              {9, 13, 341}, {10, 11, 835}, {11, 16, 316}, {12, 15, 948}, {12, 16, 40},
              {13, 16, 727} }, 1552 }
);

template <class Algorithm>
void expect_mst_can_run_twice() {
    auto G = build_graph(
        4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, false);
    auto mst = Algorithm(G);

    for (int iteration = 0; iteration < 2; ++iteration) {
        mst.run();
        int tree_weight = 0;
        mst.getForest().forEdges(
            [&](NetworKit::node, NetworKit::node, NetworKit::edgeweight weight) {
                tree_weight += weight;
            });
        EXPECT_EQ(mst.getForest().numberOfEdges(), 3);
        EXPECT_EQ(tree_weight, 20);
        mst.check();
    }
}

void expect_approximate_mst_can_run_twice() {
    auto G = build_graph(
        4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, false);
    auto mst = Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree(G);
    constexpr unsigned int max_weight = 15;
    constexpr float eps = 0.45;

    for (int iteration = 0; iteration < 2; ++iteration) {
        mst.run(max_weight, eps);
        EXPECT_GE(mst.getTreeWeight(), 20 * (1 - eps));
        EXPECT_LE(mst.getTreeWeight(), 20 * (1 + eps));
    }
}

template <class Algorithm>
void expect_mst_rejects_directed_graph() {
    auto G = build_graph(2, {{0, 1, 1}}, true);
    auto mst = Algorithm(G);

    EXPECT_THROW(mst.run(), std::runtime_error);
}

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
            [&](NetworKit::node, NetworKit::node, NetworKit::edgeweight w) {
                tree_weight += w;
            });
        EXPECT_EQ(parameters.N - 1, mst.getForest().numberOfEdges());
        EXPECT_EQ(tree_weight, parameters.treeWeight);
        mst.check();
    }

    void test_approximate_mst() {
        SpanningTreeParameters const& parameters = GetParam();
        auto G = build_graph(parameters.N, parameters.EW, false);
        auto mst = Algorithm(G);

        int max_w = 0;
        G.forEdges([&max_w](NetworKit::node, NetworKit::node,
            NetworKit::edgeweight ew, NetworKit::edgeid) {
            max_w = std::max(max_w, static_cast<int>(ew));
        });

        if (max_w == 0) {
            EXPECT_EQ(parameters.treeWeight, 0);
            return;
        }

        for (float eps : {0.005f, 0.01f, 0.025f, 0.05f, 0.15f, 0.25f, 0.35f, 0.45f}) {
            mst.run(max_w, eps);
            EXPECT_GE(mst.getTreeWeight(), parameters.treeWeight * (1 - eps));
            EXPECT_LE(mst.getTreeWeight(), parameters.treeWeight * (1 + eps));
        }
    }
};

class KruskalMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::KruskalMinimumSpanningTree> { };

TEST_P(KruskalMinimumSpanningTreeTest, test_example) {
    test_mst();
}

TEST(KruskalMinimumSpanningTreeTest, can_run_twice) {
    expect_mst_can_run_twice<Koala::KruskalMinimumSpanningTree>();
}

TEST(KruskalMinimumSpanningTreeTest, rejects_directed_graph) {
    expect_mst_rejects_directed_graph<Koala::KruskalMinimumSpanningTree>();
}

INSTANTIATE_TEST_SUITE_P(test_example_graphs, KruskalMinimumSpanningTreeTest, example_graphs);

class PrimMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::PrimMinimumSpanningTree> { };

TEST_P(PrimMinimumSpanningTreeTest, test_example) {
    test_mst();
}

TEST(PrimMinimumSpanningTreeTest, can_run_twice) {
    expect_mst_can_run_twice<Koala::PrimMinimumSpanningTree>();
}

TEST(PrimMinimumSpanningTreeTest, rejects_directed_graph) {
    expect_mst_rejects_directed_graph<Koala::PrimMinimumSpanningTree>();
}

INSTANTIATE_TEST_SUITE_P(test_example, PrimMinimumSpanningTreeTest, example_graphs);

class BoruvkaMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::BoruvkaMinimumSpanningTree> { };

TEST_P(BoruvkaMinimumSpanningTreeTest, test_example) {
    test_mst();
}

TEST(BoruvkaMinimumSpanningTreeTest, can_run_twice) {
    expect_mst_can_run_twice<Koala::BoruvkaMinimumSpanningTree>();
}

TEST(BoruvkaMinimumSpanningTreeTest, rejects_directed_graph) {
    expect_mst_rejects_directed_graph<Koala::BoruvkaMinimumSpanningTree>();
}

INSTANTIATE_TEST_SUITE_P(test_example, BoruvkaMinimumSpanningTreeTest, example_graphs);

class KargerKleinTarjanMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::KargerKleinTarjanMinimumSpanningTree> { };

TEST_P(KargerKleinTarjanMinimumSpanningTreeTest, test_example) {
    test_mst();
}

TEST(KargerKleinTarjanMinimumSpanningTreeTest, can_run_twice) {
    expect_mst_can_run_twice<Koala::KargerKleinTarjanMinimumSpanningTree>();
}

TEST(KargerKleinTarjanMinimumSpanningTreeTest, rejects_directed_graph) {
    expect_mst_rejects_directed_graph<Koala::KargerKleinTarjanMinimumSpanningTree>();
}

INSTANTIATE_TEST_SUITE_P(test_example, KargerKleinTarjanMinimumSpanningTreeTest, example_graphs);

class ChazelleRubinfeldTrevisanMinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree> { };

TEST_P(ChazelleRubinfeldTrevisanMinimumSpanningTreeTest, test_example) {
    test_approximate_mst();
}

INSTANTIATE_TEST_SUITE_P(test_example, ChazelleRubinfeldTrevisanMinimumSpanningTreeTest,
    example_graphs);

TEST(ChazelleRubinfeldTrevisanMinimumSpanningTreeTest, can_run_twice) {
    expect_approximate_mst_can_run_twice();
}

TEST(ChazelleRubinfeldTrevisanMinimumSpanningTreeTest, rejects_directed_graph) {
    auto G = build_graph(2, {{0, 1, 1}}, true);
    auto mst = Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree(G);

    EXPECT_THROW(mst.run(1), std::runtime_error);
}

class Chazelle2000MinimumSpanningTreeTest
    : public MinimumSpanningTreeTest<Koala::Chazelle2000MinimumSpanningTree> { };

TEST_P(Chazelle2000MinimumSpanningTreeTest, test_example) {
    test_mst();
}

TEST(Chazelle2000MinimumSpanningTreeTest, can_run_twice) {
    expect_mst_can_run_twice<Koala::Chazelle2000MinimumSpanningTree>();
}

TEST(Chazelle2000MinimumSpanningTreeTest, rejects_directed_graph) {
    expect_mst_rejects_directed_graph<Koala::Chazelle2000MinimumSpanningTree>();
}

INSTANTIATE_TEST_CASE_P(test_example, Chazelle2000MinimumSpanningTreeTest, example_graphs);
