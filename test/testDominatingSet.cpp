#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <list>

#include <networkit/graph/AdjListGraph.hpp>

#include <dominating_set/BakerKOuterplanarGraphDominatingSet.hpp>
#include <dominating_set/BakerPlanarGraphDominatingSet.hpp>
#include <dominating_set/ExactDominatingSet.hpp>
#include <set_cover/BranchAndReduceSetCover.hpp>

#include "test/helpers.hpp"

struct DominatingSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int minimumDominatingSetSize;
};

namespace {

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

void verifyBakerDominatingSet(
        NetworKit::Graph graph, std::size_t expected_size) {
    Koala::BakerKOuterplanarGraphDominatingSet algorithm(graph);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(algorithm.getDominatingSet().size(), expected_size)
        << "solution="
        << ::testing::PrintToString(algorithm.getDominatingSet())
        << " levels="
        << ::testing::PrintToString(algorithm.getBakerForest().levels);
    EXPECT_TRUE(algorithm.getBakerForest().hasTwoKBoundaryBound());
}

}  // namespace

class GrandoniTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class FominGrandoniKratschTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class RooijBodlaenderTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class FominKratschWoegingerTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class SchiermeyerTest
    : public testing::TestWithParam<DominatingSetParameters> {};

class BakerGraph6ErrorsDominatingSetTest
    : public testing::TestWithParam<DominatingSetParameters> {};

auto parameter_set = DominatingSetParameters{
    32,
    {{0, 1}, {0, 2}, {0, 4}, {0, 8}, {0, 16}, {1, 2}, {1, 3}, {1, 5}, {1, 9}, {1, 17},
        {2, 3}, {2, 4}, {2, 6}, {2, 10}, {2, 18}, {3, 4}, {3, 5}, {3, 7}, {3, 11}, {3, 19},
        {4, 5}, {4, 6}, {4, 8}, {4, 12}, {4, 20}, {5, 6}, {5, 7}, {5, 9}, {5, 13}, {5, 21},
        {6, 7}, {6, 8}, {6, 10}, {6, 14}, {6, 22}, {7, 8}, {7, 9}, {7, 11}, {7, 15}, {7, 23},
        {8, 9}, {8, 10}, {8, 12}, {8, 16}, {8, 24}, {9, 10}, {9, 11}, {9, 13}, {9, 17}, {9, 25},
        {10, 11}, {10, 12}, {10, 14}, {10, 18}, {10, 26}, {11, 12}, {11, 13}, {11, 15},
        {11, 19}, {11, 27}, {12, 13}, {12, 14}, {12, 16}, {12, 20}, {12, 28}, {13, 14},
        {13, 15}, {13, 17}, {13, 21}, {13, 29}, {14, 15}, {14, 16}, {14, 18}, {14, 22},
        {14, 30}, {15, 16}, {15, 17}, {15, 19}, {15, 23}, {15, 31}, {16, 17}, {16, 18},
        {16, 20}, {16, 24}, {17, 18}, {17, 19}, {17, 21}, {17, 25}, {18, 19}, {18, 20},
        {18, 22}, {18, 26}, {19, 20}, {19, 21}, {19, 23}, {19, 27}, {20, 21}, {20, 22},
        {20, 24}, {20, 28}, {21, 22}, {21, 23}, {21, 25}, {21, 29}, {22, 23}, {22, 24},
        {22, 26}, {22, 30}, {23, 24}, {23, 25}, {23, 27}, {23, 31}, {24, 25}, {24, 26},
        {24, 28}, {25, 26}, {25, 27}, {25, 29}, {26, 27}, {26, 28}, {26, 30}, {27, 28},
        {27, 29}, {27, 31}, {28, 29}, {28, 30}, {29, 30}, {29, 31}, {30, 31}},
    5};

TEST_P(GrandoniTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, GrandoniTest, testing::Values(parameter_set));

TEST_P(FominGrandoniKratschTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, FominGrandoniKratschTest, testing::Values(parameter_set));

TEST_P(RooijBodlaenderTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, RooijBodlaenderTest, testing::Values(parameter_set));

TEST_P(FominKratschWoegingerTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::FominKratschWoegingerDominatingSet(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, FominKratschWoegingerTest, testing::Values(parameter_set));

TEST_P(SchiermeyerTest, test) {
    DominatingSetParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::SchiermeyerDominatingSet(G);
    algorithm.run();
    algorithm.check();
    EXPECT_EQ(parameters.minimumDominatingSetSize, algorithm.getDominatingSet().size());
}

INSTANTIATE_TEST_SUITE_P(test_example, SchiermeyerTest, testing::Values(parameter_set));

TEST_P(BakerGraph6ErrorsDominatingSetTest, SolvesRegressionGraph) {
    const auto &parameters = GetParam();
    const auto graph = build_graph(parameters.N, parameters.E, false);
    verifyBakerDominatingSet(
        graph, parameters.minimumDominatingSetSize);
}

INSTANTIATE_TEST_SUITE_P(
    BakerErrors, BakerGraph6ErrorsDominatingSetTest, testing::Values(
    // G@P\~{
    DominatingSetParameters{
        8,
        {{2, 3}, {1, 4}, {1, 5}, {3, 5}, {4, 5},
         {0, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6},
         {0, 7}, {1, 7}, {2, 7}, {3, 7}, {4, 7}, {5, 7}, {6, 7}},
        1},
    // H??ZVIV
    DominatingSetParameters{
        9,
        {{3, 5}, {4, 5}, {1, 6}, {2, 6}, {4, 6},
         {0, 7}, {1, 7}, {2, 7}, {5, 7},
         {0, 8}, {3, 8}, {5, 8}, {6, 8}, {7, 8}},
        2},
    // I???FE]^_
    DominatingSetParameters{
        10,
        {{0, 7}, {1, 7}, {2, 7}, {6, 7},
         {0, 8}, {3, 8}, {4, 8}, {5, 8}, {6, 8},
         {1, 9}, {2, 9}, {3, 9}, {4, 9}, {5, 9}, {6, 9}},
        2}));

TEST(BakerKOuterplanarGraphDominatingSetTest, TriangleChains) {
    constexpr std::array<std::size_t, 7> expected_sizes = {
        1, 2, 3, 4, 4, 5, 6
    };
    for (NetworKit::count number_of_triangles = 1;
         number_of_triangles <= expected_sizes.size();
         ++number_of_triangles) {
        SCOPED_TRACE(number_of_triangles);
        verifyBakerDominatingSet(
            makeTriangleChain(number_of_triangles),
            expected_sizes[number_of_triangles - 1]);
    }
}

TEST(BakerKOuterplanarGraphDominatingSetTest, SquareChains) {
    constexpr std::array<std::size_t, 7> expected_sizes = {
        2, 2, 3, 4, 5, 6, 7
    };
    for (NetworKit::count number_of_squares = 1;
         number_of_squares <= expected_sizes.size();
         ++number_of_squares) {
        SCOPED_TRACE(number_of_squares);
        verifyBakerDominatingSet(
            makeSquareChain(number_of_squares),
            expected_sizes[number_of_squares - 1]);
    }
}

TEST(BakerPlanarGraphDominatingSetTest, TriangleAndSquareChains) {
    constexpr NetworKit::count parameter = 1;
    constexpr std::array<std::size_t, 3> triangle_optima = {1, 2, 3};
    constexpr std::array<std::size_t, 3> square_optima = {2, 2, 3};
    for (NetworKit::count number_of_faces = 1;
         number_of_faces <= triangle_optima.size(); ++number_of_faces) {
        for (const auto &[graph, optimum] : {
                 std::pair{
                     makeTriangleChain(number_of_faces),
                     triangle_optima[number_of_faces - 1]},
                 std::pair{
                     makeSquareChain(number_of_faces),
                     square_optima[number_of_faces - 1]}}) {
            SCOPED_TRACE(number_of_faces);
            NetworKit::Graph mutable_graph = graph;
            Koala::BakerPlanarGraphDominatingSet algorithm(
                mutable_graph, 1.0);
            algorithm.run();
            algorithm.check();

            const auto result_size = algorithm.getDominatingSet().size();
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

TEST(BakerPlanarGraphDominatingSetTest, UsesInheritedOuterFace) {
    auto graph = makeOuterFaceRegressionGraph();
    Koala::BakerPlanarGraphDominatingSet algorithm(graph, 1.0);
    algorithm.run();
    algorithm.check();

    const auto &statistics = algorithm.getBakerApproximationStatistics();
    EXPECT_EQ(statistics.parameter, 1);
    EXPECT_EQ(statistics.promisedOuterplanarity, 2);
    EXPECT_LE(
        statistics.maximumSubproblemOuterplanarity,
        statistics.promisedOuterplanarity);
}
