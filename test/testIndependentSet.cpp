#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <list>
#include <set>
#include <vector>

#include "independent_set/BakerKOuterplanarGraphIndependentSet.hpp"
#include "independent_set/BakerPlanarGraphIndependentSet.hpp"
#include "independent_set/ChordalIndependentSet.hpp"
#include "independent_set/CographIndependentSet.hpp"
#include "independent_set/IndependentSet.hpp"
#include "recognition/ChordalGraphRecognition.hpp"
#include "recognition/CographRecognition.hpp"

#include "test/helpers.hpp"

struct IndependentSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

template <typename Algorithm>
class SimpleGraphs : public testing::Test {
 protected:
    void verify(IndependentSetParameters& parameters) {
        NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
        auto algorithm = Algorithm(G);
        algorithm.run();

        auto independentSet = algorithm.getIndependentSet();
        std::set<NetworKit::node> iset(independentSet.begin(), independentSet.end());
        for (const auto &[u, v] : parameters.E) {
            EXPECT_FALSE(iset.contains(u) && iset.contains(v));
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
        {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4},
        {5, 6}, {5, 7}, {5, 8}, {5, 9}, {6, 7}, {6, 8}, {6, 9}, {7, 8}, {7, 9}, {8, 9},
    },
    2};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, TestingGraph) {
    IndependentSetParameters parameters =
    {6, {
        {1, 2}, {0, 3}, {1, 3}, {2, 3}, {0, 4}, {2, 4}, {0, 5}, {1, 5},
    },
    3};
    this->verify(parameters);
}
REGISTER_TYPED_TEST_CASE_P(
    SimpleGraphs, WheelGraphW_8, UtilityGraphK_3_3, PetersenGraph, FruchtGraph,
    TwoK5, TestingGraph);

using Algorithms = testing::Types<
    Koala::BruteForceIndependentSet,
    Koala::Mis1IndependentSet,
    Koala::Mis2IndependentSet,
    Koala::Mis3IndependentSet,
    Koala::Mis4IndependentSet,
    Koala::Mis5IndependentSet,
    Koala::MeasureAndConquerIndependentSet>;

INSTANTIATE_TYPED_TEST_CASE_P(IndependentSet, SimpleGraphs, Algorithms);

class ChordalIndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };

TEST_P(ChordalIndependentSetTest, test) {
    IndependentSetParameters const &parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto recognition = Koala::MaximumCardinalitySearchChordalGraphRecognition(G);
    recognition.run();
    EXPECT_TRUE(recognition.isChordal());
    auto algorithm = Koala::ChordalIndependentSet(G, recognition.getPEO());
    algorithm.run();
    auto iset = algorithm.getIndependentSet();
    std::set<NetworKit::node> nodes(iset.begin(), iset.end());
    for (const auto &[u, v] : parameters.E)
        EXPECT_FALSE(nodes.contains(u) && nodes.contains(v));
    EXPECT_EQ(static_cast<int>(iset.size()), parameters.expectedSetSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, ChordalIndependentSetTest, testing::Values(
    IndependentSetParameters{0, {}, 0},
    IndependentSetParameters{1, {}, 1},
    IndependentSetParameters{4, {{0, 1}, {0, 2}, {0, 3}}, 3},
    IndependentSetParameters{4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}}, 2},
    IndependentSetParameters{
        5, {{0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 1},
    IndependentSetParameters{5, {{0, 1}, {1, 2}, {2, 0}, {3, 4}}, 2},
    IndependentSetParameters{
        6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
            {1, 4}, {1, 5}, {2, 4}}, 2}
));

TEST(CographIndependentSetTest, mapsCotreeLeavesToGraphNodes) {
    IndependentSetParameters parameters = {
        6, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {5, 1}, {5, 2}, {5, 3}, {5, 4}}, 4};
    auto G = build_graph(parameters.N, parameters.E, false);
    Koala::HabibPaulCographRecognition recognition(G);
    recognition.run();
    ASSERT_TRUE(recognition.isCograph());

    Koala::CographIndependentSet algorithm(G, recognition.cotree);
    algorithm.run();
    auto independent_set = algorithm.getIndependentSet();
    std::set<NetworKit::node> nodes(independent_set.begin(), independent_set.end());
    for (auto node : independent_set) {
        EXPECT_TRUE(G.hasNode(node));
    }
    for (const auto &[u, v] : parameters.E) {
        EXPECT_FALSE(nodes.contains(u) && nodes.contains(v));
    }
    EXPECT_EQ(independent_set.size(), parameters.expectedSetSize);
}

namespace {

std::size_t bakerPower(std::size_t base, std::size_t exponent) {
    std::size_t result = 1;
    for (std::size_t i = 0; i < exponent; ++i) {
        result *= base;
    }
    return result;
}

NetworKit::Graph makeTriangleChain(NetworKit::count number_of_triangles) {
    NetworKit::Graph graph(3 * number_of_triangles);
    for (NetworKit::count i = 0; i < number_of_triangles; ++i) {
        graph.addEdge(3 * i, 3 * i + 1);
        graph.addEdge(3 * i + 1, 3 * i + 2);
        graph.addEdge(3 * i + 2, 3 * i);
    }
    for (NetworKit::node j = 0; j + 3 < 3 * number_of_triangles; ++j) {
        graph.addEdge(j, j + 3);
    }
    return graph;
}

NetworKit::Graph makeSquareChain(NetworKit::count number_of_squares) {
    NetworKit::Graph graph(4 * number_of_squares);
    for (NetworKit::count i = 0; i < number_of_squares; ++i) {
        graph.addEdge(4 * i, 4 * i + 1);
        graph.addEdge(4 * i + 1, 4 * i + 2);
        graph.addEdge(4 * i + 2, 4 * i + 3);
        graph.addEdge(4 * i + 3, 4 * i);
    }
    for (NetworKit::node j = 0; j + 4 < 4 * number_of_squares; ++j) {
        graph.addEdge(j, j + 4);
    }
    return graph;
}

NetworKit::Graph makeOuterFaceRegressionGraph() {
    // G??KzW
    NetworKit::Graph graph(8);
    constexpr std::array<std::pair<NetworKit::node, NetworKit::node>, 9>
        edges = {{
            {4, 5},
            {0, 6}, {3, 6}, {4, 6}, {5, 6},
            {1, 7}, {2, 7}, {4, 7}, {5, 7}
        }};
    for (const auto &[u, v] : edges) {
        graph.addEdge(u, v);
    }
    return graph;
}

NetworKit::Graph makeWheelEight() {
    NetworKit::Graph graph(8);
    for (NetworKit::node vertex = 0; vertex < 7; ++vertex) {
        graph.addEdge(vertex, (vertex + 1) % 7);
        graph.addEdge(vertex, 7);
    }
    return graph;
}

NetworKit::Graph makeConcentricCycles(NetworKit::count cycle_size) {
    NetworKit::Graph graph(2 * cycle_size);
    for (NetworKit::node vertex = 0; vertex < cycle_size; ++vertex) {
        const auto next = (vertex + 1) % cycle_size;
        graph.addEdge(vertex, next);
        graph.addEdge(cycle_size + vertex, cycle_size + next);
        graph.addEdge(vertex, cycle_size + vertex);
    }
    return graph;
}

void verifyBakerForest(const Koala::BakerForest &forest) {
    EXPECT_TRUE(forest.hasTwoKBoundaryBound());
    EXPECT_EQ(forest.width(), 2 * forest.outerplanarity());

    std::size_t number_of_tree_vertices = 0;
    std::vector<bool> is_root(forest.trees.size(), false);
    for (const auto root : forest.roots) {
        ASSERT_LT(root, forest.trees.size());
        is_root[root] = true;
    }

    for (std::size_t tree_index = 0;
         tree_index < forest.trees.size(); ++tree_index) {
        const auto &tree = forest.trees[tree_index];
        number_of_tree_vertices += tree.vertices.size();
        ASSERT_FALSE(tree.vertices.empty());
        ASSERT_LT(tree.root, tree.vertices.size());
        EXPECT_FALSE(tree.vertices[tree.root].parent.has_value());
        EXPECT_EQ(is_root[tree_index], !tree.enclosingTree.has_value());
        EXPECT_EQ(
            tree.enclosingTree.has_value(),
            tree.enclosingVertex.has_value());

        if (tree.enclosingTree.has_value()) {
            ASSERT_LT(*tree.enclosingTree, forest.trees.size());
            const auto &parent_tree = forest.trees[*tree.enclosingTree];
            ASSERT_LT(*tree.enclosingVertex, parent_tree.vertices.size());
            ASSERT_TRUE(
                parent_tree.vertices[*tree.enclosingVertex]
                    .enclosedTree.has_value());
            EXPECT_EQ(
                *parent_tree.vertices[*tree.enclosingVertex].enclosedTree,
                tree_index);
        }

        for (std::size_t vertex_index = 0;
             vertex_index < tree.vertices.size(); ++vertex_index) {
            const auto &vertex = tree.vertices[vertex_index];
            EXPECT_EQ(vertex.boundary.left.size(), tree.level);
            EXPECT_EQ(vertex.boundary.right.size(), tree.level);
            EXPECT_LE(vertex.boundary.size(), 2 * forest.outerplanarity());
            EXPECT_LE(
                vertex.leftBoundaryIndex, vertex.rightBoundaryIndex);
            ASSERT_FALSE(vertex.boundary.left.empty());
            ASSERT_FALSE(vertex.boundary.right.empty());
            EXPECT_EQ(vertex.boundary.left.front(), vertex.labelLeft);
            EXPECT_EQ(vertex.boundary.right.front(), vertex.labelRight);

            for (std::size_t i = 0; i < tree.level; ++i) {
                const auto expected_level = tree.level - i;
                ASSERT_LT(vertex.boundary.left[i], forest.levels.size());
                ASSERT_LT(vertex.boundary.right[i], forest.levels.size());
                EXPECT_EQ(
                    forest.levels[vertex.boundary.left[i]],
                    expected_level);
                EXPECT_EQ(
                    forest.levels[vertex.boundary.right[i]],
                    expected_level);
            }

            if (vertex.parent.has_value()) {
                ASSERT_LT(*vertex.parent, tree.vertices.size());
                const auto &siblings =
                    tree.vertices[*vertex.parent].children;
                EXPECT_NE(
                    std::find(
                        siblings.begin(), siblings.end(), vertex_index),
                    siblings.end());
            } else {
                EXPECT_EQ(vertex_index, tree.root);
            }

            for (std::size_t child = 0;
                 child < vertex.children.size(); ++child) {
                const auto child_index = vertex.children[child];
                ASSERT_LT(child_index, tree.vertices.size());
                ASSERT_TRUE(tree.vertices[child_index].parent.has_value());
                EXPECT_EQ(
                    *tree.vertices[child_index].parent, vertex_index);
                if (child > 0) {
                    const auto previous = vertex.children[child - 1];
                    EXPECT_EQ(
                        tree.vertices[previous].boundary.right,
                        tree.vertices[child_index].boundary.left);
                }
            }
            if (vertex.type == Koala::BakerFaceTreeVertexType::FACE) {
                ASSERT_FALSE(vertex.children.empty());
                EXPECT_EQ(
                    vertex.labelLeft,
                    tree.vertices[vertex.children.front()].labelLeft);
                EXPECT_EQ(
                    vertex.labelRight,
                    tree.vertices[vertex.children.back()].labelRight);
            } else {
                EXPECT_TRUE(vertex.children.empty());
                EXPECT_LE(
                    vertex.leftBoundaryIndex,
                    vertex.middleBoundaryIndex);
                EXPECT_LE(
                    vertex.middleBoundaryIndex,
                    vertex.rightBoundaryIndex);
            }
            if (vertex.enclosedTree.has_value()) {
                ASSERT_LT(*vertex.enclosedTree, forest.trees.size());
                EXPECT_EQ(
                    forest.trees[*vertex.enclosedTree].enclosingTree,
                    tree_index);
                EXPECT_EQ(
                    forest.trees[*vertex.enclosedTree].enclosingVertex,
                    vertex_index);
            }
        }

        for (std::size_t leaf = 0; leaf < tree.leaves.size(); ++leaf) {
            ASSERT_LT(tree.leaves[leaf], tree.vertices.size());
            const auto &current = tree.vertices[tree.leaves[leaf]];
            EXPECT_EQ(
                current.type,
                Koala::BakerFaceTreeVertexType::EXTERIOR_EDGE);
            if (leaf > 0) {
                const auto &previous =
                    tree.vertices[tree.leaves[leaf - 1]];
                EXPECT_EQ(previous.labelRight, current.labelLeft);
            }
        }
    }

    const auto &statistics = forest.statistics;
    EXPECT_EQ(statistics.recursiveCalls, number_of_tree_vertices);
    EXPECT_LE(statistics.adjustCalls, number_of_tree_vertices);
    EXPECT_LE(statistics.mergeCalls, number_of_tree_vertices);
    EXPECT_LE(statistics.contractCalls, number_of_tree_vertices);
    EXPECT_LE(statistics.createCalls, number_of_tree_vertices);
    EXPECT_LE(statistics.extendCalls, number_of_tree_vertices);

    const std::size_t table_bound =
        bakerPower(2, 2 * forest.outerplanarity());
    const std::size_t merge_bound =
        bakerPower(2, 3 * forest.outerplanarity());
    const std::size_t operation_calls =
        statistics.adjustCalls + statistics.mergeCalls
        + statistics.contractCalls + statistics.createCalls
        + statistics.extendCalls;
    EXPECT_LE(statistics.maximumTableEntries, table_bound);
    EXPECT_LE(statistics.maximumMergeTransitions, merge_bound);
    EXPECT_LE(
        statistics.tableEntries,
        (number_of_tree_vertices + operation_calls) * table_bound);
    EXPECT_LE(
        statistics.mergeTransitions,
        statistics.mergeCalls * merge_bound);
}

}  // namespace

void verifyBakerIndependentSet(
        const NetworKit::Graph &graph, NetworKit::count expected_size) {
    Koala::BakerKOuterplanarGraphIndependentSet algorithm(graph);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(algorithm.getIndependentSet().size(), expected_size);
    verifyBakerForest(algorithm.getBakerForest());
}

class BakerGraph6ErrorsIndependentSetTest
    : public testing::TestWithParam<IndependentSetParameters> { };

TEST_P(BakerGraph6ErrorsIndependentSetTest, SolvesRegressionGraph) {
    const auto &parameters = GetParam();
    const auto graph = build_graph(parameters.N, parameters.E, false);
    verifyBakerIndependentSet(graph, parameters.expectedSetSize);
}

INSTANTIATE_TEST_SUITE_P(
    BakerErrors, BakerGraph6ErrorsIndependentSetTest, testing::Values(
    // E@lw
    IndependentSetParameters{
        6,
        {{2, 3}, {0, 4}, {2, 4}, {3, 4},
         {1, 5}, {2, 5}, {3, 5}, {4, 5}},
        3},
    // F?C}o
    IndependentSetParameters{
        7,
        {{3, 4}, {2, 5}, {3, 5}, {4, 5},
         {0, 6}, {1, 6}, {3, 6}, {4, 6}},
        4},
    // G??B~{
    IndependentSetParameters{
        8,
        {{1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6},
         {0, 7}, {1, 7}, {2, 7}, {3, 7}, {4, 7}, {5, 7}, {6, 7}},
        6},
    // H???@~~
    IndependentSetParameters{
        9,
        {{2, 7}, {3, 7}, {4, 7}, {5, 7}, {6, 7},
         {0, 8}, {1, 8}, {2, 8}, {3, 8}, {4, 8}, {5, 8}, {6, 8}, {7, 8}},
        7},
    // I?????^~w
    IndependentSetParameters{
        10,
        {{3, 8}, {4, 8}, {5, 8}, {6, 8}, {7, 8},
         {0, 9}, {1, 9}, {2, 9}, {3, 9}, {4, 9}, {5, 9}, {6, 9}, {7, 9}, {8, 9}},
        8},
    // FDZ^w
    IndependentSetParameters{
        7,
        {{0, 3}, {2, 3}, {1, 4}, {2, 4}, {0, 5}, {1, 5}, {3, 5},
         {4, 5}, {0, 6}, {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6}},
        3},
    // G??ZfS
    IndependentSetParameters{
        8,
        {{3, 5}, {4, 5}, {1, 6}, {2, 6}, {3, 6},
         {0, 7}, {1, 7}, {2, 7}, {4, 7}, {6, 7}},
        5},
    // H???Xrt
    IndependentSetParameters{
        9,
        {{4, 6}, {5, 6}, {2, 7}, {3, 7}, {4, 7},
         {0, 8}, {1, 8}, {2, 8}, {3, 8}, {5, 8}, {7, 8}},
        6},
    // I????E|^g
    IndependentSetParameters{
        10,
        {{6, 7}, {0, 8}, {2, 8}, {3, 8}, {4, 8}, {5, 8}, {7, 8},
         {1, 9}, {2, 9}, {3, 9}, {4, 9}, {5, 9}, {6, 9}, {8, 9}},
        7}));

TEST(BakerKOuterplanarGraphIndependentSetTest, EmptyAndIsolatedVertices) {
    verifyBakerIndependentSet(NetworKit::Graph(0), 0);

    NetworKit::Graph isolated_vertices(5);
    Koala::BakerKOuterplanarGraphIndependentSet algorithm(isolated_vertices);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(algorithm.getIndependentSet().size(), 5);
    EXPECT_EQ(algorithm.getBakerForest().outerplanarity(), 1);
    EXPECT_EQ(algorithm.getBakerForest().roots.size(), 5);
}

TEST(BakerKOuterplanarGraphIndependentSetTest, BridgesAndFakeEdges) {
    NetworKit::Graph path(3);
    path.addEdge(0, 1);
    path.addEdge(1, 2);

    Koala::BakerKOuterplanarGraphIndependentSet algorithm(path);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(algorithm.getIndependentSet().size(), 2);

    std::size_t bridge_copy_count = 0;
    for (const auto &edge : algorithm.getBakerForest().fakeEdges) {
        if (edge.type == Koala::BakerFakeEdgeType::BRIDGE_COPY) {
            ++bridge_copy_count;
            EXPECT_TRUE(path.hasEdge(edge.u, edge.v));
        }
    }
    EXPECT_EQ(bridge_copy_count, 2);
    verifyBakerForest(algorithm.getBakerForest());
}

TEST(BakerKOuterplanarGraphIndependentSetTest, ManyBridgesRemainLinearSize) {
    constexpr NetworKit::count number_of_vertices = 256;
    NetworKit::Graph path(number_of_vertices);
    for (NetworKit::node vertex = 1;
         vertex < number_of_vertices; ++vertex) {
        path.addEdge(vertex - 1, vertex);
    }

    Koala::BakerKOuterplanarGraphIndependentSet algorithm(path);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(
        algorithm.getIndependentSet().size(),
        (number_of_vertices + 1) / 2);
    EXPECT_EQ(
        std::count_if(
            algorithm.getBakerForest().fakeEdges.begin(),
            algorithm.getBakerForest().fakeEdges.end(),
            [](const Koala::BakerFakeEdge &edge) {
                return edge.type
                    == Koala::BakerFakeEdgeType::BRIDGE_COPY;
            }),
        number_of_vertices - 1);
    verifyBakerForest(algorithm.getBakerForest());
}

TEST(BakerKOuterplanarGraphIndependentSetTest, SelfLoopIsAnOriginalConstraint) {
    NetworKit::Graph graph(2);
    graph.addEdge(0, 0);
    graph.addEdge(0, 1);
    verifyBakerIndependentSet(graph, 1);
}

TEST(BakerKOuterplanarGraphIndependentSetTest, PreservesSparseNodeIdentifiers) {
    NetworKit::Graph graph(6);
    graph.removeNode(1);
    graph.removeNode(3);
    graph.removeNode(4);
    graph.addEdge(0, 2);
    graph.addEdge(2, 5);
    verifyBakerIndependentSet(graph, 2);
}

TEST(BakerKOuterplanarGraphIndependentSetTest, ChainsOfTriangles) {
    for (NetworKit::count number_of_triangles = 1;
         number_of_triangles <= 7; ++number_of_triangles) {
        SCOPED_TRACE(number_of_triangles);
        verifyBakerIndependentSet(
            makeTriangleChain(number_of_triangles),
            number_of_triangles);
    }
}

TEST(BakerKOuterplanarGraphIndependentSetTest, ChainsOfSquares) {
    for (NetworKit::count number_of_squares = 1;
         number_of_squares <= 7; ++number_of_squares) {
        SCOPED_TRACE(number_of_squares);
        verifyBakerIndependentSet(
            makeSquareChain(number_of_squares), 2 * number_of_squares);
    }
}

TEST(
    BakerKOuterplanarGraphIndependentSetTest,
    MaterializesOrderedFaceTreesAndNamedOperations) {
    const auto graph = makeTriangleChain(7);
    Koala::BakerKOuterplanarGraphIndependentSet algorithm(graph);
    algorithm.run();
    algorithm.check();
    verifyBakerForest(algorithm.getBakerForest());

    bool has_nested_tree = false;
    bool has_nontrivial_lb_rb = false;
    for (const auto &tree : algorithm.getBakerForest().trees) {
        has_nested_tree = has_nested_tree || tree.enclosingTree.has_value();
        for (const auto &vertex : tree.vertices) {
            has_nontrivial_lb_rb =
                has_nontrivial_lb_rb
                || vertex.leftBoundaryIndex != vertex.rightBoundaryIndex;
        }
    }
    EXPECT_TRUE(has_nested_tree);
    EXPECT_TRUE(has_nontrivial_lb_rb);

    const auto &statistics = algorithm.getBakerForest().statistics;
    EXPECT_GT(statistics.adjustCalls, 0);
    EXPECT_GT(statistics.mergeCalls, 0);
    EXPECT_GT(statistics.contractCalls, 0);
    EXPECT_GT(statistics.createCalls, 0);
    EXPECT_GT(statistics.extendCalls, 0);
}

TEST(
    BakerKOuterplanarGraphIndependentSetTest,
    NestedSingletonComponentUsesExplicitTree) {
    const auto graph = makeWheelEight();
    Koala::BakerKOuterplanarGraphIndependentSet algorithm(graph);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(algorithm.getIndependentSet().size(), 3);
    verifyBakerForest(algorithm.getBakerForest());

    bool found_nested_singleton = false;
    for (const auto &tree : algorithm.getBakerForest().trees) {
        found_nested_singleton =
            found_nested_singleton
            || (tree.level == 2
                && tree.vertices.size() == 1
                && tree.vertices.front().type
                    == Koala::BakerFaceTreeVertexType::SINGLETON
                && tree.enclosingTree.has_value());
    }
    EXPECT_TRUE(found_nested_singleton);
}

TEST(
    BakerKOuterplanarGraphIndependentSetTest,
    LongAnnulusUsesMaterializedSliceMaps) {
    constexpr NetworKit::count cycle_size = 32;
    const auto graph = makeConcentricCycles(cycle_size);
    Koala::BakerKOuterplanarGraphIndependentSet algorithm(graph);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(algorithm.getIndependentSet().size(), cycle_size);
    EXPECT_EQ(algorithm.getBakerForest().outerplanarity(), 2);
    verifyBakerForest(algorithm.getBakerForest());
}

TEST(
    BakerKOuterplanarGraphIndependentSetTest,
    MatchesBruteForceOnWheelSubgraphs) {
    const std::vector<std::pair<NetworKit::node, NetworKit::node>> edges{
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 0},
        {7, 0}, {7, 1}, {7, 2}, {7, 3}, {7, 4}, {7, 5}, {7, 6}};
    for (std::size_t sample = 0; sample < 64; ++sample) {
        SCOPED_TRACE(sample);
        NetworKit::Graph graph(8);
        const std::size_t edge_mask =
            sample * 0x9e3779b97f4a7c15ULL;
        for (std::size_t edge = 0; edge < edges.size(); ++edge) {
            if ((edge_mask >> edge) % 2 == 1) {
                graph.addEdge(edges[edge].first, edges[edge].second);
            }
        }

        Koala::BruteForceIndependentSet brute_force(graph);
        brute_force.run();
        Koala::BakerKOuterplanarGraphIndependentSet baker(graph);
        baker.run();
        baker.check();
        EXPECT_EQ(
            baker.getIndependentSet().size(),
            brute_force.getIndependentSet().size());
        EXPECT_TRUE(baker.getBakerForest().hasTwoKBoundaryBound());
    }
}

TEST(BakerKOuterplanarGraphIndependentSetTest, DisconnectedPlanarGraph) {
    NetworKit::Graph graph(8);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 0);
    graph.addEdge(3, 4);
    graph.addEdge(4, 5);
    graph.addEdge(5, 6);
    graph.addEdge(6, 3);
    verifyBakerIndependentSet(graph, 4);
}

TEST(BakerKOuterplanarGraphIndependentSetTest, RejectsNonplanarGraph) {
    NetworKit::Graph graph(6);
    for (NetworKit::node first = 0; first < 3; ++first) {
        for (NetworKit::node second = 3; second < 6; ++second) {
            graph.addEdge(first, second);
        }
    }
    Koala::BakerKOuterplanarGraphIndependentSet algorithm(graph);
    EXPECT_THROW(algorithm.run(), std::invalid_argument);
}

TEST(BakerPlanarGraphIndependentSetTest, TriangleAndSquareChains) {
    constexpr NetworKit::count parameter = 2;
    for (NetworKit::count number_of_faces = 1;
         number_of_faces <= 4; ++number_of_faces) {
        for (const auto &[graph, optimum] : {
                 std::pair{makeTriangleChain(number_of_faces), number_of_faces},
                 std::pair{makeSquareChain(number_of_faces), 2 * number_of_faces}}) {
            SCOPED_TRACE(number_of_faces);
            Koala::BakerPlanarGraphIndependentSet algorithm(graph, 0.5);
            algorithm.run();
            algorithm.check();

            const auto result_size = algorithm.getIndependentSet().size();
            EXPECT_GE(result_size * (parameter + 1), optimum * parameter);
            const auto &statistics =
                algorithm.getBakerApproximationStatistics();
            EXPECT_EQ(statistics.parameter, parameter);
            EXPECT_EQ(statistics.promisedOuterplanarity, parameter);
            EXPECT_LE(
                statistics.maximumSubproblemOuterplanarity,
                statistics.promisedOuterplanarity);
            EXPECT_EQ(statistics.residueCandidates, parameter + 1);
        }
    }
}

TEST(BakerPlanarGraphIndependentSetTest, EmptyAndInvalidInput) {
    NetworKit::Graph empty_graph(0);
    Koala::BakerPlanarGraphIndependentSet empty_algorithm(empty_graph, 0.5);
    empty_algorithm.run();
    empty_algorithm.check();
    EXPECT_TRUE(empty_algorithm.getIndependentSet().empty());

    EXPECT_THROW(
        Koala::BakerPlanarGraphIndependentSet(empty_graph, 0.0),
        std::invalid_argument);

    NetworKit::Graph nonplanar_graph(6);
    for (NetworKit::node first = 0; first < 3; ++first) {
        for (NetworKit::node second = 3; second < 6; ++second) {
            nonplanar_graph.addEdge(first, second);
        }
    }
    Koala::BakerPlanarGraphIndependentSet nonplanar_algorithm(
        nonplanar_graph, 0.5);
    EXPECT_THROW(nonplanar_algorithm.run(), std::invalid_argument);
}

TEST(BakerPlanarGraphIndependentSetTest, TranslatesSparseNodeIdentifiers) {
    NetworKit::Graph graph(8);
    graph.addEdge(0, 3);
    graph.addEdge(3, 7);
    for (const NetworKit::node removed : {1, 2, 4, 5, 6}) {
        graph.removeNode(removed);
    }

    Koala::BakerPlanarGraphIndependentSet algorithm(graph, 0.5);
    algorithm.run();
    algorithm.check();
    EXPECT_GE(algorithm.getIndependentSet().size(), 2);
    for (const auto vertex : algorithm.getIndependentSet()) {
        EXPECT_TRUE(graph.hasNode(vertex));
    }
}

TEST(BakerPlanarGraphIndependentSetTest, UsesInheritedOuterFace) {
    const auto graph = makeOuterFaceRegressionGraph();
    Koala::BakerPlanarGraphIndependentSet algorithm(graph, 1.0);
    algorithm.run();
    algorithm.check();

    const auto &statistics = algorithm.getBakerApproximationStatistics();
    EXPECT_EQ(statistics.parameter, 1);
    EXPECT_EQ(statistics.promisedOuterplanarity, 1);
    EXPECT_LE(
        statistics.maximumSubproblemOuterplanarity,
        statistics.promisedOuterplanarity);
}
