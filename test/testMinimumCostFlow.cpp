// TODO: undirected graph, antisymmetric costs
#include <gtest/gtest.h>

#include <list>

#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>

#include "helpers.hpp"

struct MinimumCostCirculationParameters {
    int N;
    // from, to, capacity, cost
    std::list<std::tuple<int, int, int, int>> EW;
    int minCost;
};

class SuccessiveApproxMCCTest
    : public testing::TestWithParam<MinimumCostCirculationParameters> { };

TEST_P(SuccessiveApproxMCCTest, test) {
    MinimumCostCirculationParameters const& parameters = GetParam();
    std::list<std::pair<int,int>> edges;
    Koala::edge_map<Koala::MCFEdgeParams> ep;

    for(auto [v,w,cap, cost] : parameters.EW){
        edges.push_back({v,w});
        ep[{v,w}].capacity += cap;
        ep[{v,w}].cost += cost;
        ep[{w,v}].cost -= cost;
    }
    NetworKit::Graph G = build_graph(parameters.N, edges, false);
    auto algorithm = Koala::SuccessiveApproxMCC(G, ep);
    algorithm.run();
    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, SuccessiveApproxMCCTest, testing::Values(
        MinimumCostCirculationParameters{
            5, {{0, 1, 10, 2}, {1,3,4,3},{3,1,-2,0},{0, 2, 15, 1}, {2, 3, 3, -3}, {3, 0, 10, 0}} , 
        }
        
));
