#include <gtest/gtest.h>

#include <cmath>

#include <networkit/graph/AdjListGraph.hpp>

#include "flow/electrical_flow/FlowNetwork.hpp"

constexpr double tolerance = 1e-7;

void expect_conservation(
        const Koala::FlowNetwork &network, NetworKit::node source, NetworKit::node target,
        double value) {
    for (NetworKit::node u = 0; u < network.graph.numberOfNodes(); ++u) {
        double balance = 0;
        for (NetworKit::node v = 0; v < network.graph.numberOfNodes(); ++v) {
            EXPECT_NEAR(network.flow[u][v], -network.flow[v][u], tolerance);
            balance += network.flow[u][v];
        }
        EXPECT_NEAR(balance, u == source ? -value : u == target ? value : 0, tolerance);
    }
}

TEST(FlowNetworkTest, computes_residual_capacities) {
    NetworKit::Graph graph(2, true, true);
    graph.addEdge(0, 1, 5);
    Koala::FlowNetwork network(graph);
    network.flow[0][1] = 2;
    network.flow[1][0] = -2;

    EXPECT_DOUBLE_EQ(network.upperCapacity(0, 1), 3);
    EXPECT_DOUBLE_EQ(network.lowerCapacity(0, 1), 7);
}

TEST(FlowNetworkTest, pushes_value_through_multiple_paths) {
    NetworKit::Graph graph(4, true, false);
    graph.addEdge(0, 1, 2);
    graph.addEdge(1, 3, 2);
    graph.addEdge(0, 2, 3);
    graph.addEdge(2, 3, 3);
    Koala::FlowNetwork network(graph);

    EXPECT_TRUE(network.pushValue(0, 3, 5));

    expect_conservation(network, 0, 3, 5);
    for (auto [u, v] : graph.edgeRange()) {
        EXPECT_LE(std::abs(network.flow[u][v]), graph.weight(u, v));
    }
}

TEST(FlowNetworkTest, reports_missing_residual_path) {
    NetworKit::Graph graph(4, true, false);
    graph.addEdge(0, 1, 2);
    graph.addEdge(2, 3, 2);
    Koala::FlowNetwork network(graph);

    EXPECT_FALSE(network.pushValue(0, 3, 1));
}

TEST(FlowNetworkTest, rounds_fractional_cycle_without_changing_balance) {
    NetworKit::Graph graph(3, true, false);
    graph.addEdge(0, 1, 2);
    graph.addEdge(1, 2, 2);
    graph.addEdge(2, 0, 2);
    Koala::FlowNetwork network(graph);
    network.flow[0][1] = 0.5, network.flow[1][0] = -0.5;
    network.flow[1][2] = 0.5, network.flow[2][1] = -0.5;
    network.flow[2][0] = 0.5, network.flow[0][2] = -0.5;

    network.roundFlow();

    expect_conservation(network, 0, 0, 0);
    for (auto [u, v] : graph.edgeRange()) {
        EXPECT_NEAR(network.flow[u][v], std::round(network.flow[u][v]), tolerance);
    }
}

TEST(FlowNetworkTest, rounds_near_integer_flow_without_adding_a_unit) {
    NetworKit::Graph graph(2, true, false);
    graph.addEdge(0, 1, 1);
    Koala::FlowNetwork network(graph);
    network.flow[0][1] = -5e-9;
    network.flow[1][0] = 5e-9;

    network.roundFlow();

    EXPECT_DOUBLE_EQ(network.flow[0][1], 0);
    EXPECT_DOUBLE_EQ(network.flow[1][0], 0);
}
