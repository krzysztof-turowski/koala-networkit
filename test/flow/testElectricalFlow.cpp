#include <gtest/gtest.h>

#include <vector>

#include <networkit/graph/AdjListGraph.hpp>

#include "flow/electrical_flow/ElectricalNetwork.hpp"

constexpr double tolerance = 1e-7;

TEST(ElectricalNetworkTest, routes_demand_and_preserves_antisymmetry) {
    NetworKit::Graph graph(3, true, false);
    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 1);
    std::vector<double> demand{-3, 0, 3};
    std::vector<std::vector<double>> resistance(3, std::vector<double>(3, 0));
    resistance[0][1] = resistance[1][0] = 2;
    resistance[1][2] = resistance[2][1] = 4;
    Koala::ElectricalNetwork network(graph, demand);

    network.compute(resistance);

    for (NetworKit::node u = 0; u < graph.numberOfNodes(); ++u) {
        double balance = 0;
        for (NetworKit::node v = 0; v < graph.numberOfNodes(); ++v) {
            EXPECT_NEAR(network.flow[u][v], -network.flow[v][u], tolerance);
            balance += network.flow[u][v];
        }
        EXPECT_NEAR(balance, demand[u], tolerance);
    }
}
