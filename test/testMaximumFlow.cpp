#include <gtest/gtest.h>

#include <cmath>
#include <list>
#include <random>
#include <stdexcept>

#include "flow/MaximumFlow.hpp"
#include "flow/GoldbergTarjanPushRelabelMaximumFlow.hpp"
#include "flow/BoykovKolmogorovFlow.hpp"
#include "flow/electrical_flow/ElectricalFlow.hpp"
#include "flow/electrical_flow/FlowNetwork.hpp"
#include "flow/KingRaoTarjanMaximumFlow.hpp"
#include "flow/MalhotraKumarMaheshwariFlow.hpp"
#include "flow/maximum_flow/KrtEdgeDesignator.hpp"

#include "test/helpers.hpp"

struct MaximumFlowParameters {
    int N;
    std::list<std::tuple<int, int, int>> EW;
    int s, t;
    int flowSize;
};

template <class Algorithm>
void expect_maximum_flow_can_run_twice() {
    NetworKit::Graph G = build_graph(
        4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, true);
    auto algorithm = Algorithm(G, 0, 3);

    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), 15);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), 15);
}

class ElectricalFlowTest : public testing::TestWithParam<MaximumFlowParameters> { };

NetworKit::Graph build_electrical_graph(const MaximumFlowParameters &parameters) {
    return build_graph(parameters.N, parameters.EW, true);
}

void expect_electrical_flow(
        const Koala::ElectricalFlow &algorithm, NetworKit::node source, NetworKit::node target) {
    const auto &G = algorithm.getGraph();
    const auto &flow = algorithm.getFlow();
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        if (G.isDirected()) {
            EXPECT_GE(flow[u][v], -G.weight(u, v));
        } else {
            EXPECT_LE(std::abs(flow[u][v]), G.weight(u, v));
        }
        EXPECT_NEAR(flow[u][v], std::round(flow[u][v]), 1e-7);
        EXPECT_NEAR(flow[u][v], -flow[v][u], 1e-7);
    });
    G.forNodes([&](NetworKit::node u) {
        double balance = 0;
        for (NetworKit::node v = 0; v < G.numberOfNodes(); ++v) {
            balance += flow[u][v];
        }
        EXPECT_NEAR(
            balance, u == source ? -algorithm.getFlowSize()
                                 : u == target ? algorithm.getFlowSize() : 0, 1e-7);
    });
}

TEST_P(ElectricalFlowTest, computes_maximum_flow) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_electrical_graph(parameters);
    auto algorithm = Koala::ElectricalFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
    expect_electrical_flow(algorithm, parameters.s, parameters.t);
}

INSTANTIATE_TEST_SUITE_P(
    execution_paths, ElectricalFlowTest, testing::Values(
        MaximumFlowParameters{2, {}, 0, 1, 0},
        MaximumFlowParameters{2, {{0, 1, 7}}, 0, 1, 7},
        MaximumFlowParameters{4, {{0, 1, 5}, {1, 3, 5}, {0, 2, 7}, {2, 3, 7}}, 0, 3, 12},
        MaximumFlowParameters{
            4, {{0, 1, 10}, {1, 2, 3}, {2, 3, 10}, {0, 3, 1}}, 0, 3, 4},
        MaximumFlowParameters{
            4, {{0, 1, 3}, {0, 2, 5}, {1, 2, 2}, {2, 1, 3}, {1, 3, 7}, {2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{4, {{0, 1, 5}, {2, 3, 5}}, 0, 3, 0}));

TEST(ElectricalFlowTest, can_run_twice) {
    NetworKit::Graph G = build_electrical_graph(
        {4, {{0, 1, 5}, {1, 3, 5}, {0, 2, 7}, {2, 3, 7}}, 0, 3, 12});
    auto algorithm = Koala::ElectricalFlow(G, 0, 3);

    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), 12);
    expect_electrical_flow(algorithm, 0, 3);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), 12);
    expect_electrical_flow(algorithm, 0, 3);
}

TEST(ElectricalFlowTest, computes_undirected_maximum_flow) {
    NetworKit::Graph G = build_graph(
        4, {{0, 1, 5}, {1, 3, 5}, {0, 2, 7}, {2, 3, 7}}, false);
    auto algorithm = Koala::ElectricalFlow(G, 0, 3);

    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), 12);
    expect_electrical_flow(algorithm, 0, 3);
}

TEST(FlowNetworkTest, rounds_near_integer_flow_without_adding_a_unit) {
    NetworKit::Graph G(2, true, false);
    G.addEdge(0, 1, 1);
    Koala::FlowNetwork network(G);
    network.flow[0][1] = -5e-9;
    network.flow[1][0] = 5e-9;

    network.roundFlow();

    EXPECT_DOUBLE_EQ(network.flow[0][1], 0);
    EXPECT_DOUBLE_EQ(network.flow[1][0], 0);
}

class KingRaoTarjanMaximumFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(KingRaoTarjanMaximumFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::KingRaoTarjanMaximumFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, KingRaoTarjanMaximumFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4, {{0, 1, 3}, {0, 2, 5}, {1, 2, 2}, {2, 1, 3}, {1, 3, 7}, {2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {{0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                 {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                 {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                 {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}}, 0, 21, 0}
));

TEST(KingRaoTarjanMaximumFlowTest, test_prim_designator_nodes) {
    constexpr NetworKit::count paths = 512;
    constexpr NetworKit::node source = paths, target = paths + 1;
    NetworKit::Graph G(paths + 2, true, true);
    for (NetworKit::node v = 0; v < paths; ++v) {
        G.addEdge(source, v, 1);
        G.addEdge(v, target, 1);
    }

    auto algorithm = Koala::KingRaoTarjanMaximumFlow(
        G, source, target, {.L = 512, .T = 7});
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), paths);
}

TEST(KingRaoTarjanMaximumFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::KingRaoTarjanMaximumFlow>();
}

TEST(KingRaoTarjanMaximumFlowTest, matches_push_relabel_on_small_random_graphs) {
    std::mt19937 generator(123456);
    for (NetworKit::count n = 2; n <= 8; ++n) {
        for (int iteration = 0; iteration < 20; ++iteration) {
            NetworKit::Graph G(n, true, true);
            for (NetworKit::node u = 0; u < n; ++u) {
                for (NetworKit::node v = 0; v < n; ++v) {
                    if (u != v && generator() % 4 == 0) {
                        G.addEdge(u, v, 1 + generator() % 20);
                    }
                }
            }

            auto expected = Koala::GoldbergTarjanPushRelabelMaximumFlow(G, 0, n - 1);
            expected.run();
            auto actual = Koala::KingRaoTarjanMaximumFlow(G, 0, n - 1);
            actual.run();
            EXPECT_EQ(actual.getFlowSize(), expected.getFlowSize())
                << "n = " << n << ", iteration = " << iteration;
        }
    }
}

TEST(KRTEdgeDesignatorTest, redesignates_after_designated_edge_removal) {
    NetworKit::Graph G(3, true, true);
    G.addEdge(0, 1, 1);
    G.addEdge(0, 2, 1);
    Koala::KRTEdgeDesignator designator;
    designator.initialize(G, {.L = 512, .T = 7});

    auto first = designator.current_edge(0, 1);
    ASSERT_NE(first, NetworKit::none);
    designator.response_adversary(0, 1, first, 0);

    auto second = designator.current_edge(0, 1);
    EXPECT_NE(second, NetworKit::none);
    EXPECT_NE(second, first);
    designator.response_adversary(0, 1, second, 0);
    EXPECT_EQ(designator.current_edge(0, 1), NetworKit::none);
}

TEST(KRTEdgeDesignatorTest, redesignates_after_leaving_U_prim) {
    constexpr NetworKit::count neighbors = 512;
    NetworKit::Graph G(neighbors + 1, true, true);
    for (NetworKit::node v = 1; v <= neighbors; ++v) {
        G.addEdge(0, v, 1);
    }
    Koala::KRTEdgeDesignator designator;
    designator.initialize(G, {.L = 512, .T = 7});

    auto first = designator.current_edge(0, 1);
    ASSERT_NE(first, NetworKit::none);
    designator.response_adversary(0, 1, first, 0);

    auto second = designator.current_edge(0, 1);
    EXPECT_NE(second, NetworKit::none);
    EXPECT_NE(second, first);
}

TEST(KRTEdgeDesignatorTest, rejects_invalid_parameters) {
    NetworKit::Graph G(2, true, true);
    G.addEdge(0, 1, 1);
    Koala::KRTEdgeDesignator designator;

    EXPECT_THROW(designator.initialize(G, {.L = 1}), std::invalid_argument);
}

class GoldbergTarjanPushRelabelMaximumFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(GoldbergTarjanPushRelabelMaximumFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm =
        Koala::GoldbergTarjanPushRelabelMaximumFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, GoldbergTarjanPushRelabelMaximumFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4, {{0, 1, 3}, {0, 2, 5}, {1, 2, 2}, {2, 1, 3}, {1, 3, 7}, {2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {{0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                 {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                 {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                 {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}}, 0, 21, 0}
));

TEST(GoldbergTarjanPushRelabelMaximumFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::GoldbergTarjanPushRelabelMaximumFlow>();
}

class MKMFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(MKMFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::MalhotraKumarMaheshwariFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, MKMFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4, {{0, 1, 3}, {0, 2, 5}, {1, 2, 2}, {2, 1, 3}, {1, 3, 7}, {2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {{0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                 {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                 {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                 {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}}, 0, 21, 0}
));

TEST(MKMFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::MalhotraKumarMaheshwariFlow>();
}

class BKFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(BKFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::BoykovKolmogorovFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, BKFlowTest, testing::Values(
        MaximumFlowParameters{
            4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
        MaximumFlowParameters{
            6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
                {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
        MaximumFlowParameters{
            4, {{0, 1, 3}, {0, 2, 5}, {1, 2, 2}, {2, 1, 3}, {1, 3, 7}, {2, 3, 1}}, 0, 3, 7},
        MaximumFlowParameters{
            22, {{0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
                 {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
                 {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
                 {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}}, 0, 21, 0}
));

TEST(BKFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::BoykovKolmogorovFlow>();
}
