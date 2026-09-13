#include <gtest/gtest.h>

#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include <networkit/distance/APSP.hpp>
#include <networkit/graph/Graph.hpp>

#include "shortest_path/Aingworth.hpp"
#include "shortest_path/Seidel.hpp"
#include "shortest_path/UnweightedShoshanZwick.hpp"
#include "shortest_path/Spira.hpp"

namespace {

constexpr int INF = std::numeric_limits<int>::max();

NetworKit::Graph unweighted(int n, const std::vector<std::pair<int, int>>& edges) {
    NetworKit::Graph G(n, false, false);
    for (const auto& [u, v] : edges) {
        G.addEdge(u, v);
    }
    return G;
}

NetworKit::Graph weighted(int n, const std::vector<std::tuple<int, int, double>>& edges) {
    NetworKit::Graph G(n, true, false);
    for (const auto& [u, v, w] : edges) {
        G.addEdge(u, v, w);
    }
    return G;
}

NetworKit::Graph grid(int rows, int cols) {
    NetworKit::Graph G(rows * cols, false, false);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const int node = r * cols + c;
            if (c + 1 < cols) {
                G.addEdge(node, node + 1);
            }
            if (r + 1 < rows) {
                G.addEdge(node, node + cols);
            }
        }
    }
    return G;
}

template <typename Algorithm>
void expectExactUnweighted(NetworKit::Graph& G) {
    NetworKit::APSP reference(G);
    reference.run();
    Algorithm algorithm(G);
    algorithm.run();
    for (const auto u : G.nodeRange()) {
        for (const auto v : G.nodeRange()) {
            EXPECT_EQ(algorithm.getDistance(u, v),
                      static_cast<int>(reference.getDistance(u, v)))
                << "d(" << u << ", " << v << ")";
        }
    }
}

}  // namespace

struct APSPGraph {
    std::string name;
    NetworKit::Graph (*build)();
};

class APSPTest : public testing::TestWithParam<APSPGraph> {};

TEST_P(APSPTest, SeidelMatchesReference) {
    NetworKit::Graph G = GetParam().build();
    expectExactUnweighted<Koala::SeidelAPSP>(G);
}

TEST_P(APSPTest, UnweightedShoshanZwickMatchesReference) {
    NetworKit::Graph G = GetParam().build();
    expectExactUnweighted<Koala::UnweightedShoshanZwickAPSP>(G);
}

// Aingworth is an additive +2 approximation: true <= estimate <= true + 2.
TEST_P(APSPTest, AingworthWithinAdditiveTwo) {
    NetworKit::Graph G = GetParam().build();
    NetworKit::APSP reference(G);
    reference.run();
    Koala::AingworthAPSP algorithm(G);
    algorithm.run();
    for (const auto u : G.nodeRange()) {
        for (const auto v : G.nodeRange()) {
            const int exact = static_cast<int>(reference.getDistance(u, v));
            const int estimate = algorithm.getDistance(u, v);
            EXPECT_GE(estimate, exact) << "d(" << u << ", " << v << ")";
            EXPECT_LE(estimate, exact + 2) << "d(" << u << ", " << v << ")";
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Default, APSPTest,
    testing::Values(
        APSPGraph{"path", [] { return unweighted(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}}); }},
        APSPGraph{"cycle",
            [] { return unweighted(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}); }},
        APSPGraph{"star", [] { return unweighted(6, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}}); }},
        APSPGraph{"grid3x3", [] { return grid(3, 3); }},
        APSPGraph{"general",
            [] { return unweighted(7, {{0, 1}, {0, 2}, {1, 3}, {2, 3}, {3, 4}, {4, 5}, {4, 6}}); }},
        APSPGraph{"grid5x5", [] { return grid(5, 5); }}),
    [](const testing::TestParamInfo<APSPGraph>& info) { return info.param.name; });

// Spira handles weighted graphs: check it against NetworKit::APSP (Dijkstra).
TEST(APSPWeightedTest, SpiraMatchesReference) {
    NetworKit::Graph G = weighted(5,
        {{0, 1, 2.0}, {1, 2, 1.0}, {0, 3, 5.0}, {3, 2, 1.0}, {2, 4, 3.0}});
    NetworKit::APSP reference(G);
    reference.run();
    Koala::SpiraAPSP algorithm(G);
    algorithm.run();
    for (const auto u : G.nodeRange()) {
        for (const auto v : G.nodeRange()) {
            EXPECT_DOUBLE_EQ(algorithm.getDistance(u, v), reference.getDistance(u, v))
                << "d(" << u << ", " << v << ")";
        }
    }
}

// Trivial edge cases: a single node and a single edge.
TEST(APSPTest, TrivialGraphs) {
    NetworKit::Graph single(1, false, false);
    Koala::SeidelAPSP a(single);
    a.run();
    EXPECT_EQ(a.getDistance(0, 0), 0);
    EXPECT_EQ(a.getDiameter(), 0);

    NetworKit::Graph edge = unweighted(2, {{0, 1}});
    Koala::UnweightedShoshanZwickAPSP b(edge);
    b.run();
    EXPECT_EQ(b.getDistance(0, 1), 1);
    EXPECT_EQ(b.getDistance(0, 0), 0);
}

// Non-consecutive node ids (a deleted node) must not break the matrix algorithms.
TEST(APSPTest, NonConsecutiveNodeIds) {
    NetworKit::Graph G(5, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 3);
    G.addEdge(3, 4);
    G.removeNode(2);
    Koala::SeidelAPSP algorithm(G);
    algorithm.run();
    EXPECT_EQ(algorithm.getDistance(0, 4), 3);
    EXPECT_EQ(algorithm.getDiameter(), 3);
}

// Input-validation preconditions surface as exceptions.
TEST(APSPTest, InputValidationThrows) {
    NetworKit::Graph directed(3, false, true);
    EXPECT_THROW(Koala::SeidelAPSP{directed}, std::invalid_argument);

    NetworKit::Graph weightedGraph(3, true, false);
    EXPECT_THROW(Koala::AingworthAPSP{weightedGraph}, std::invalid_argument);

    NetworKit::Graph unweightedGraph(3, false, false);
    EXPECT_THROW(Koala::SpiraAPSP{unweightedGraph}, std::invalid_argument);

    // The matrix algorithms require a connected graph.
    NetworKit::Graph disconnected = unweighted(4, {{0, 1}, {2, 3}});
    EXPECT_THROW(Koala::SeidelAPSP{disconnected}, std::invalid_argument);
    EXPECT_THROW(Koala::UnweightedShoshanZwickAPSP{disconnected}, std::invalid_argument);
}

// Aingworth and Spira leave unreachable pairs at infinity on disconnected input.
TEST(APSPDisconnectedTest, AingworthAndSpiraReturnInfinity) {
    NetworKit::Graph G = unweighted(4, {{0, 1}, {2, 3}});
    Koala::AingworthAPSP aingworth(G);
    aingworth.run();
    EXPECT_EQ(aingworth.getDistance(0, 2), INF);
    EXPECT_EQ(aingworth.getDistance(0, 1), 1);

    NetworKit::Graph W = weighted(4, {{0, 1, 2.0}, {2, 3, 5.0}});
    Koala::SpiraAPSP spira(W);
    spira.run();
    EXPECT_EQ(spira.getDistance(0, 2), std::numeric_limits<NetworKit::edgeweight>::max());
    EXPECT_DOUBLE_EQ(spira.getDistance(0, 1), 2.0);
}
