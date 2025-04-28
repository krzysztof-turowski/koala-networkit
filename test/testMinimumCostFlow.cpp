// TODO: undirected graph, antisymmetric costs
#include <gtest/gtest.h>

#include <list>

#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>

#include "helpers.hpp"

struct MinCostCirculationParams {
    int N;
    // from, to, capacity, cost
    std::list<std::tuple<int, int, int, int>> EW;
    int minCost;
};

struct MinCostFlowParams {
    int N;
    // from, to, capacity, cost
    std::list<std::tuple<int, int, int, int>> EW;
    std::unordered_map<int,int> excess;
    int minCost;
};

struct MinCostFlowParamsST {
    int N;
    // from, to, capacity, cost
    std::list<std::tuple<int, int, int, int>> EW;
    int s,t;
    int flow;
};

class SuccessiveApproxMCCTest
    : public testing::TestWithParam<MinCostCirculationParams> { };

class SuccessiveApproxMCFlowTest
    : public testing::TestWithParam<MinCostFlowParams> { };
 

TEST_P(SuccessiveApproxMCCTest, test) {
    MinCostCirculationParams const& parameters = GetParam();
    std::list<std::pair<int,int>> edges;
    Koala::edge_map<Koala::MCFEdgeParams> ep;

    for(auto [v,w,cap, cost] : parameters.EW){
        edges.push_back({v,w});
        ep[{v,w}].capacity += cap;
        ep[{v,w}].cost += cost;
        ep[{w,v}].cost -= cost;
    }
    NetworKit::Graph G = build_graph(parameters.N, edges, false);
    G.removeMultiEdges();
    auto algorithm = Koala::SuccessiveApproxMCC(G, ep);
    algorithm.run();
    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

TEST_P(SuccessiveApproxMCFlowTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    std::list<std::pair<int,int>> edges;
    Koala::edge_map<Koala::MCFEdgeParams> ep;
    Koala::node_map<int> excess;
    for(auto [a,b] : parameters.excess)
        excess[a]=b;

    for(auto [v,w,cap, cost] : parameters.EW){
        edges.push_back({v,w});
        ep[{v,w}].capacity += cap;
        ep[{v,w}].cost += cost;
        ep[{w,v}].cost -= cost;
    }
    NetworKit::Graph G = build_graph(parameters.N, edges, false);
    G.removeMultiEdges();
    auto algorithm = Koala::SuccessiveApproxMCC(G, ep, excess);
    algorithm.run();
    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}


INSTANTIATE_TEST_SUITE_P(
    test_example, SuccessiveApproxMCCTest, testing::Values(
        MinCostCirculationParams{
            5, {{0, 1, 10, 2}, {1,3,4,3}, {3,1,-2,0}, {0, 2, 15, 1}, {2, 3, 3, -3}, {3, 0, 10, 0}}, 4
        }
        
));

INSTANTIATE_TEST_SUITE_P(
    test_flow, SuccessiveApproxMCFlowTest, testing::Values(
        MinCostFlowParams{
            5, {{0, 1, 10, 2}, {1,3,4,3}, {3,1,-2,0}, {0, 2, 15, 1}, {2, 3, 3, -3}},{{0,3},{3,-3}}, 8
        }
));
