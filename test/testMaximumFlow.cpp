#include <gtest/gtest.h>

#include <cmath>
#include <list>
#include <random>

#include "flow/MaximumFlow.hpp"
#include "flow/GoldbergTarjanPushRelabelMaximumFlow.hpp"
#include "flow/BoykovKolmogorovMaximumFlow.hpp"
#include "flow/ElectricalFlow.hpp"
#include "flow/KingRaoTarjanMaximumFlow.hpp"
#include "flow/MalhotraKumarMaheshwariFlow.hpp"

#include "test/helpers.hpp"

struct MaximumFlowParameters {
    int N;
    std::list<std::tuple<int, int, int>> EW;
    int s, t;
    int flowSize;
};

auto example_graphs = testing::Values(
    MaximumFlowParameters{2, {}, 0, 1, 0},
    MaximumFlowParameters{2, {{0, 1, 7}}, 0, 1, 7},
    MaximumFlowParameters{4, {{0, 1, 5}, {1, 3, 5}, {0, 2, 7}, {2, 3, 7}}, 0, 3, 12},
    MaximumFlowParameters{
        4, {{0, 1, 10}, {1, 2, 3}, {2, 3, 10}, {0, 3, 1}}, 0, 3, 4},
    MaximumFlowParameters{
        4, {{0, 1, 3}, {0, 2, 5}, {1, 2, 2}, {2, 1, 3}, {1, 3, 7}, {2, 3, 1}}, 0, 3, 7},
    MaximumFlowParameters{4, {{0, 1, 5}, {2, 3, 5}}, 0, 3, 0},
    MaximumFlowParameters{
        4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, 0, 3, 15},
    MaximumFlowParameters{
        6, {{0, 1, 16}, {0, 2, 13}, {1, 3, 12}, {2, 1, 4}, {2, 4, 14}, {3, 2, 9},
            {4, 3, 7}, {3, 5, 20}, {4, 5, 4}}, 0, 5, 23},
    MaximumFlowParameters{
        22, {{0, 1, 10}, {0, 2, 10}, {0, 3, 10}, {0, 4, 10}, {0, 5, 10},
             {0, 6, 10}, {0, 7, 10}, {0, 8, 10}, {0, 9, 10}, {0, 10, 10},
             {11, 21, 10}, {12, 21, 10}, {13, 21, 10}, {14, 21, 10}, {15, 21, 10},
             {16, 21, 10}, {17, 21, 10}, {18, 21, 10}, {19, 21, 10}, {20, 21, 10}},
        0, 21, 0});

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

template <class Algorithm>
void expect_maximum_flow_converts_undirected_graph() {
    NetworKit::Graph G = build_graph(
        4, {{0, 1, 10}, {0, 2, 5}, {1, 2, 15}, {1, 3, 5}, {2, 3, 10}}, false);
    auto algorithm = Algorithm(G, 0, 3);

    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), 15);
}

class ElectricalFlowTest : public testing::TestWithParam<MaximumFlowParameters> { };

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
        if (u == source) {
            EXPECT_NEAR(balance, -algorithm.getFlowSize(), 1e-7);
        } else if (u == target) {
            EXPECT_NEAR(balance, algorithm.getFlowSize(), 1e-7);
        } else {
            EXPECT_NEAR(balance, 0, 1e-7);
        }
    });
}

TEST_P(ElectricalFlowTest, computes_maximum_flow) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::ElectricalFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
    expect_electrical_flow(algorithm, parameters.s, parameters.t);
}

INSTANTIATE_TEST_SUITE_P(test_example, ElectricalFlowTest, example_graphs);

TEST(ElectricalFlowTest, can_run_twice) {
    NetworKit::Graph G = build_graph(
        4, {{0, 1, 5}, {1, 3, 5}, {0, 2, 7}, {2, 3, 7}}, true);
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

class KingRaoTarjanMaximumFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(KingRaoTarjanMaximumFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::KingRaoTarjanMaximumFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(test_example, KingRaoTarjanMaximumFlowTest, example_graphs);

TEST(KingRaoTarjanMaximumFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::KingRaoTarjanMaximumFlow>();
}

TEST(KingRaoTarjanMaximumFlowTest, converts_undirected_graph) {
    expect_maximum_flow_converts_undirected_graph<Koala::KingRaoTarjanMaximumFlow>();
}

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

INSTANTIATE_TEST_SUITE_P(test_example, GoldbergTarjanPushRelabelMaximumFlowTest, example_graphs);

TEST(GoldbergTarjanPushRelabelMaximumFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::GoldbergTarjanPushRelabelMaximumFlow>();
}

TEST(GoldbergTarjanPushRelabelMaximumFlowTest, converts_undirected_graph) {
    expect_maximum_flow_converts_undirected_graph<
        Koala::GoldbergTarjanPushRelabelMaximumFlow>();
}

class MalhotraKumarMaheshwariFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(MalhotraKumarMaheshwariFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::MalhotraKumarMaheshwariFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(test_example, MalhotraKumarMaheshwariFlowTest, example_graphs);

TEST(MalhotraKumarMaheshwariFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::MalhotraKumarMaheshwariFlow>();
}

TEST(MalhotraKumarMaheshwariFlowTest, converts_undirected_graph) {
    expect_maximum_flow_converts_undirected_graph<Koala::MalhotraKumarMaheshwariFlow>();
}

class BoykovKolmogorovMaximumFlowTest
    : public testing::TestWithParam<MaximumFlowParameters> { };

TEST_P(BoykovKolmogorovMaximumFlowTest, test) {
    MaximumFlowParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.EW, true);
    auto algorithm = Koala::BoykovKolmogorovMaximumFlow(G, parameters.s, parameters.t);
    algorithm.run();
    EXPECT_EQ(algorithm.getFlowSize(), parameters.flowSize);
}

INSTANTIATE_TEST_SUITE_P(test_example, BoykovKolmogorovMaximumFlowTest, example_graphs);

TEST(BoykovKolmogorovMaximumFlowTest, can_run_twice) {
    expect_maximum_flow_can_run_twice<Koala::BoykovKolmogorovMaximumFlow>();
}

TEST(BoykovKolmogorovMaximumFlowTest, converts_undirected_graph) {
    expect_maximum_flow_converts_undirected_graph<Koala::BoykovKolmogorovMaximumFlow>();
}
