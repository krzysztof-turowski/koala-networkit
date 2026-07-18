#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <list>
#include <utility>
#include <vector>

#include <networkit/graph/AdjListGraph.hpp>

#include "dominating_set/BakerDominatingSet.hpp"
#include "dominating_set/BakerKOuterplanarGraphDominatingSet.hpp"
#include "dominating_set/ExactDominatingSet.hpp"
#include "independent_set/IndependentSet.hpp"
#include "vertex_cover/BakerKOuterplanarGraphVertexCover.hpp"
#include "vertex_cover/BakerPlanarGraphVertexCover.hpp"

namespace Koala {
namespace {

struct VertexCoverParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

std::size_t checkedPower(std::size_t base, std::size_t exponent) {
    std::size_t result = 1;
    for (std::size_t i = 0; i < exponent; ++i) {
        EXPECT_LE(result, std::numeric_limits<std::size_t>::max() / base);
        result *= base;
    }
    return result;
}

template <typename Algorithm>
void expectBakerBounds(const Algorithm &algorithm, std::size_t states) {
    const BakerForest &forest = algorithm.getBakerForest();
    const std::size_t k = forest.outerplanarity();
    EXPECT_TRUE(forest.hasTwoKBoundaryBound());
    EXPECT_LE(forest.width(), 2 * k);
    EXPECT_LE(
        forest.statistics.maximumTableEntries,
        checkedPower(states, 2 * k));
    EXPECT_LE(
        forest.statistics.maximumMergeTransitions,
        checkedPower(states, 3 * k));
}

void expectVertexCover(
        const NetworKit::Graph &graph, std::size_t expected_size) {
    BakerKOuterplanarGraphVertexCover vertex_cover(graph);
    vertex_cover.run();
    EXPECT_EQ(vertex_cover.getVertexCover().size(), expected_size);
    ::testing::Message cover_message;
    for (const NetworKit::node u : vertex_cover.getVertexCover()) {
        cover_message << u << " ";
    }
    bool valid_vertex_cover = true;
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        const bool covered = vertex_cover.getVertexCover().contains(u)
            || vertex_cover.getVertexCover().contains(v);
        EXPECT_TRUE(covered)
            << "uncovered edge " << u << " " << v
            << ", levels "
            << vertex_cover.getBakerForest().levels[u] << " "
            << vertex_cover.getBakerForest().levels[v]
            << ", cover " << cover_message;
        valid_vertex_cover = valid_vertex_cover && covered;
    });
    EXPECT_TRUE(valid_vertex_cover);
    if (valid_vertex_cover) {
        vertex_cover.check();
    }
    expectBakerBounds(vertex_cover, 2);
}

void expectOptima(
        NetworKit::Graph graph, std::size_t vertex_cover_size,
        std::size_t dominating_set_size) {
    expectVertexCover(graph, vertex_cover_size);

    BakerKOuterplanarGraphDominatingSet dominating_set(graph);
    dominating_set.run();
    EXPECT_EQ(dominating_set.getDominatingSet().size(), dominating_set_size);
    bool valid_dominating_set = true;
    graph.forNodes([&](NetworKit::node u) {
        bool dominated = dominating_set.getDominatingSet().contains(u);
        graph.forNeighborsOf(u, [&](NetworKit::node v) {
            dominated = dominated
                || dominating_set.getDominatingSet().contains(v);
        });
        EXPECT_TRUE(dominated) << "undominated vertex " << u;
        valid_dominating_set = valid_dominating_set && dominated;
    });
    EXPECT_TRUE(valid_dominating_set);
    if (valid_dominating_set) {
        dominating_set.check();
    }
    expectBakerBounds(dominating_set, 3);
}

void expectMatchesExistingExactAlgorithms(NetworKit::Graph graph) {
    BruteForceIndependentSet independent_set(graph);
    independent_set.run();
    const std::size_t expected_vertex_cover =
        graph.numberOfNodes()
        - independent_set.getIndependentSet().size();

    FominKratschWoegingerDominatingSet dominating_set(graph);
    dominating_set.run();
    const std::size_t expected_dominating_set =
        dominating_set.getDominatingSet().size();

    expectOptima(
        std::move(graph), expected_vertex_cover,
        expected_dominating_set);
}

NetworKit::Graph makePath(std::size_t size) {
    NetworKit::Graph graph(size);
    for (std::size_t i = 1; i < size; ++i) {
        graph.addEdge(i - 1, i);
    }
    return graph;
}

NetworKit::Graph makeCycle(std::size_t size) {
    NetworKit::Graph graph = makePath(size);
    if (size > 2) {
        graph.addEdge(size - 1, 0);
    }
    return graph;
}

NetworKit::Graph makeTriangleChain(NetworKit::count number_of_triangles) {
    NetworKit::Graph graph(3 * number_of_triangles);
    for (NetworKit::count i = 0; i < number_of_triangles; ++i) {
        graph.addEdge(3 * i, 3 * i + 1);
        graph.addEdge(3 * i + 1, 3 * i + 2);
        graph.addEdge(3 * i + 2, 3 * i);
    }
    for (NetworKit::node j = 0;
         j + 3 < 3 * number_of_triangles; ++j) {
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
    for (NetworKit::node j = 0;
         j + 4 < 4 * number_of_squares; ++j) {
        graph.addEdge(j, j + 4);
    }
    return graph;
}

NetworKit::Graph makeOuterFaceRegressionGraph() {
    // G?NP}{
    NetworKit::Graph graph(8);
    constexpr std::array<std::pair<NetworKit::node, NetworKit::node>, 15>
        edges = {{
            {2, 4}, {3, 4},
            {0, 5}, {1, 5}, {3, 5},
            {2, 6}, {3, 6}, {4, 6}, {5, 6},
            {0, 7}, {1, 7}, {3, 7}, {4, 7}, {5, 7}, {6, 7}
        }};
    for (const auto &[u, v] : edges) {
        graph.addEdge(u, v);
    }
    return graph;
}

class BakerGraph6ErrorsVertexCoverTest
    : public testing::TestWithParam<VertexCoverParameters> { };

TEST_P(BakerGraph6ErrorsVertexCoverTest, SolvesRegressionGraph) {
    const auto &parameters = GetParam();
    NetworKit::Graph graph(parameters.N);
    for (const auto &[u, v] : parameters.E) {
        graph.addEdge(u, v);
    }
    expectVertexCover(graph, parameters.expectedSetSize);
}

INSTANTIATE_TEST_SUITE_P(
    BakerErrors, BakerGraph6ErrorsVertexCoverTest, testing::Values(
    // F@X^w
    VertexCoverParameters{
        7,
        {{2, 3}, {1, 4}, {2, 4}, {1, 5}, {3, 5}, {4, 5},
         {0, 6}, {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6}},
        4},
    // G?C[Z[
    VertexCoverParameters{
        8,
        {{3, 4}, {3, 5}, {4, 5}, {0, 6}, {4, 6}, {5, 6},
         {1, 7}, {2, 7}, {4, 7}, {5, 7}, {6, 7}},
        4},
    // H??G[Lv
    VertexCoverParameters{
        9,
        {{4, 5}, {4, 6}, {5, 6}, {0, 7}, {5, 7}, {6, 7},
         {1, 8}, {2, 8}, {3, 8}, {5, 8}, {6, 8}, {7, 8}},
        4},
    // I???GMB]w
    VertexCoverParameters{
        10,
        {{5, 6}, {5, 7}, {6, 7}, {0, 8}, {6, 8}, {7, 8},
         {1, 9}, {2, 9}, {3, 9}, {4, 9}, {6, 9}, {7, 9}, {8, 9}},
        4},
    // F?Mzw
    VertexCoverParameters{
        7,
        {{2, 4}, {3, 4}, {0, 5}, {2, 5}, {3, 5}, {4, 5},
         {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6}},
        3},
    // G?DP^s
    VertexCoverParameters{
        8,
        {{3, 4}, {1, 5}, {3, 5}, {2, 6}, {4, 6}, {5, 6},
         {0, 7}, {1, 7}, {2, 7}, {3, 7}, {4, 7}, {6, 7}},
        4},
    // H???[|n
    VertexCoverParameters{
        9,
        {{4, 6}, {5, 6}, {0, 7}, {3, 7}, {4, 7}, {5, 7}, {6, 7},
         {1, 8}, {2, 8}, {4, 8}, {5, 8}, {6, 8}, {7, 8}},
        3},
    // I???GNJNo
    VertexCoverParameters{
        10,
        {{5, 6}, {5, 7}, {6, 7}, {0, 8}, {1, 8}, {4, 8}, {6, 8}, {7, 8},
         {2, 9}, {3, 9}, {4, 9}, {5, 9}, {6, 9}, {7, 9}},
        4},
    // GGC^b[
    VertexCoverParameters{
        8,
        {{1, 2}, {3, 4}, {3, 5}, {4, 5}, {0, 6}, {1, 6}, {2, 6},
         {3, 6}, {1, 7}, {2, 7}, {4, 7}, {5, 7}, {6, 7}},
        5},
    // H?@@O}|
    VertexCoverParameters{
        9,
        {{1, 5}, {2, 6}, {4, 6}, {3, 7}, {4, 7}, {5, 7}, {6, 7},
         {0, 8}, {2, 8}, {3, 8}, {4, 8}, {5, 8}, {7, 8}},
        4},
    // I??@?gNvg
    VertexCoverParameters{
        10,
        {{2, 6}, {3, 7}, {5, 7}, {4, 8}, {5, 8}, {6, 8}, {7, 8},
         {0, 9}, {1, 9}, {3, 9}, {4, 9}, {5, 9}, {6, 9}, {8, 9}},
        4}));

TEST(BakerKOuterplanarGraphCoverTest, DominatingAdapterHasThreeStates) {
    BakerDominatingSet problem;
    EXPECT_EQ(problem.stateCount(), 3);

    std::size_t first = 0;
    std::size_t second = 0;
    EXPECT_TRUE(problem.splitMergeState(
        BakerDominatingSet::SELECTED, first, second));
    EXPECT_EQ(first, BakerDominatingSet::SELECTED);
    EXPECT_EQ(second, BakerDominatingSet::SELECTED);
    EXPECT_TRUE(problem.splitMergeState(
        BakerDominatingSet::CERTIFIED, first, second));
    EXPECT_EQ(first, BakerDominatingSet::CERTIFIED);
    EXPECT_EQ(second, BakerDominatingSet::PENDING);
    EXPECT_TRUE(problem.splitMergeState(
        BakerDominatingSet::PENDING, first, second));
    EXPECT_EQ(first, BakerDominatingSet::PENDING);
    EXPECT_EQ(second, BakerDominatingSet::CERTIFIED);
}

TEST(BakerKOuterplanarGraphCoverTest, EmptyAndIsolatedVertices) {
    expectOptima(NetworKit::Graph(0), 0, 0);
    expectOptima(NetworKit::Graph(5), 0, 5);
}

TEST(BakerKOuterplanarGraphCoverTest, PathsCyclesAndBridges) {
    expectOptima(makePath(1), 0, 1);
    expectOptima(makePath(5), 2, 2);
    expectOptima(makeCycle(3), 2, 1);
    expectOptima(makeCycle(6), 3, 2);
}

TEST(BakerKOuterplanarGraphCoverTest, TriangleChains) {
    for (NetworKit::count number_of_triangles = 1;
         number_of_triangles <= 7; ++number_of_triangles) {
        SCOPED_TRACE(number_of_triangles);
        expectVertexCover(
            makeTriangleChain(number_of_triangles),
            2 * number_of_triangles);
    }
}

TEST(BakerKOuterplanarGraphCoverTest, SquareChains) {
    for (NetworKit::count number_of_squares = 1;
         number_of_squares <= 7; ++number_of_squares) {
        SCOPED_TRACE(number_of_squares);
        expectVertexCover(
            makeSquareChain(number_of_squares),
            2 * number_of_squares);
    }
}

TEST(BakerKOuterplanarGraphCoverTest, StarGraph) {
    NetworKit::Graph graph(6);
    for (NetworKit::node leaf = 1; leaf < 6; ++leaf) {
        graph.addEdge(0, leaf);
    }
    expectOptima(std::move(graph), 1, 1);
}

TEST(BakerKOuterplanarGraphCoverTest, DisconnectedGraph) {
    NetworKit::Graph graph(6);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 0);
    graph.addEdge(3, 4);
    expectOptima(std::move(graph), 3, 3);
}

TEST(BakerKOuterplanarGraphCoverTest, WheelUsesNestedSlices) {
    NetworKit::Graph graph(8);
    for (NetworKit::node rim = 1; rim < 8; ++rim) {
        graph.addEdge(0, rim);
        graph.addEdge(rim, rim == 7 ? 1 : rim + 1);
    }
    expectOptima(std::move(graph), 5, 1);
}

TEST(BakerKOuterplanarGraphCoverTest, SelfLoopIsOriginal) {
    NetworKit::Graph graph(1);
    graph.addEdge(0, 0);
    expectOptima(std::move(graph), 1, 1);
}

TEST(BakerKOuterplanarGraphCoverTest, PreservesSparseNodeIdentifiers) {
    NetworKit::Graph graph(8);
    graph.addEdge(0, 3);
    graph.addEdge(3, 7);
    for (const NetworKit::node removed : {1, 2, 4, 5, 6}) {
        graph.removeNode(removed);
    }
    expectOptima(std::move(graph), 1, 1);
}

TEST(BakerKOuterplanarGraphCoverTest, MatchesExistingExactAlgorithms) {
    const std::vector<std::pair<NetworKit::node, NetworKit::node>> edges = {
        {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7},
        {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 1}
    };
    constexpr std::uint64_t edge_mask = (std::uint64_t{1} << 14) - 1;
    for (std::uint64_t sample = 0; sample < 64; ++sample) {
        NetworKit::Graph graph(8);
        const std::uint64_t mask =
            (sample * 4051 + sample * sample * 17 + 23) & edge_mask;
        for (std::size_t edge = 0; edge < edges.size(); ++edge) {
            if ((mask & (std::uint64_t{1} << edge)) != 0) {
                graph.addEdge(edges[edge].first, edges[edge].second);
            }
        }

        SCOPED_TRACE(
            ::testing::Message() << "sample=" << sample << " mask=" << mask);
        expectMatchesExistingExactAlgorithms(std::move(graph));
    }
}

TEST(BakerKOuterplanarGraphCoverTest, MatchesExactAlgorithmsOnNestedGraphs) {
    const std::vector<std::pair<NetworKit::node, NetworKit::node>> edges = {
        {0, 1}, {1, 2}, {2, 0},
        {3, 0}, {3, 1}, {3, 2},
        {4, 0}, {4, 1}, {4, 3},
        {5, 1}, {5, 2}, {5, 3},
        {6, 0}, {6, 2}, {6, 3},
        {7, 0}, {7, 1}, {7, 4}
    };
    constexpr std::uint64_t edge_mask =
        (std::uint64_t{1} << 18) - 1;
    for (std::uint64_t sample = 0; sample <= 64; ++sample) {
        NetworKit::Graph graph(8);
        const std::uint64_t mask = sample == 64
            ? edge_mask
            : (sample * 7919 + sample * sample * 131 + 47) & edge_mask;
        for (std::size_t edge = 0; edge < edges.size(); ++edge) {
            if ((mask & (std::uint64_t{1} << edge)) != 0) {
                graph.addEdge(edges[edge].first, edges[edge].second);
            }
        }

        SCOPED_TRACE(
            ::testing::Message() << "sample=" << sample << " mask=" << mask);
        expectMatchesExistingExactAlgorithms(std::move(graph));
    }
}

TEST(BakerKOuterplanarGraphCoverTest, RejectsNonplanarGraph) {
    NetworKit::Graph graph(6);
    for (NetworKit::node left = 0; left < 3; ++left) {
        for (NetworKit::node right = 3; right < 6; ++right) {
            graph.addEdge(left, right);
        }
    }

    BakerKOuterplanarGraphVertexCover vertex_cover(graph);
    EXPECT_THROW(vertex_cover.run(), std::invalid_argument);
    BakerKOuterplanarGraphDominatingSet dominating_set(graph);
    EXPECT_THROW(dominating_set.run(), std::invalid_argument);
}

TEST(BakerPlanarGraphVertexCoverTest, TriangleAndSquareChains) {
    constexpr NetworKit::count parameter = 2;
    for (NetworKit::count number_of_faces = 1;
         number_of_faces <= 4; ++number_of_faces) {
        for (const auto &[graph, optimum] : {
                 std::pair{
                     makeTriangleChain(number_of_faces),
                     2 * number_of_faces},
                 std::pair{
                     makeSquareChain(number_of_faces),
                     2 * number_of_faces}}) {
            SCOPED_TRACE(number_of_faces);
            BakerPlanarGraphVertexCover algorithm(graph, 0.5);
            algorithm.run();
            algorithm.check();

            const auto result_size = algorithm.getVertexCover().size();
            EXPECT_LE(result_size * parameter, optimum * (parameter + 1));
            const auto &statistics =
                algorithm.getBakerApproximationStatistics();
            EXPECT_EQ(statistics.parameter, parameter);
            EXPECT_EQ(statistics.promisedOuterplanarity, parameter + 1);
            EXPECT_LE(
                statistics.maximumSubproblemOuterplanarity,
                statistics.promisedOuterplanarity);
            EXPECT_EQ(statistics.residueCandidates, parameter);
        }
    }
}

TEST(BakerPlanarGraphVertexCoverTest, UsesInheritedOuterFace) {
    const auto graph = makeOuterFaceRegressionGraph();
    BakerPlanarGraphVertexCover algorithm(graph, 1.0);
    algorithm.run();
    algorithm.check();

    const auto &statistics = algorithm.getBakerApproximationStatistics();
    EXPECT_EQ(statistics.parameter, 1);
    EXPECT_EQ(statistics.promisedOuterplanarity, 2);
    EXPECT_LE(
        statistics.maximumSubproblemOuterplanarity,
        statistics.promisedOuterplanarity);
}

}  // namespace
}  // namespace Koala
