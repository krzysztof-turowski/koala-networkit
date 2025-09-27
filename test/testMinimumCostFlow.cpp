// TODO: undirected graph, antisymmetric costs
#include <gtest/gtest.h>

#include <list>

#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>
#include <flow/minimum_cost_flow/EdmondsKarpMCF.hpp>
#include <flow/minimum_cost_flow/OrlinMCF.hpp>

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
 

// TEST_P(SuccessiveApproxMCCTest, test) {
//     MinCostCirculationParams const& parameters = GetParam();
//     std::list<std::pair<int,int>> edges;
//     Koala::edgeid_map<int> costs;
//     Koala::edgeid_map<std::pair<int,int>> bounds;

//     for(auto [v,w,cap, cost] : parameters.EW){
//         edges.push_back({v,w});
//         ep[{v,w}].capacity += cap;
//         ep[{v,w}].cost += cost;
//         ep[{w,v}].cost -= cost;
//     }
//     NetworKit::Graph G = build_graph(parameters.N, edges, false);
//     G.removeMultiEdges();
//     auto algorithm = Koala::SuccessiveApproxMCC(G, ep);
//     algorithm.run();
//     EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
// }

// TEST_P(SuccessiveApproxMCFlowTest, test) {
//     MinCostFlowParams const& parameters = GetParam();
//     std::list<std::pair<int,int>> edges;
//     Koala::edge_map<Koala::MCFEdgeParams> ep;
//     Koala::node_map<int> excess;
//     for(auto [a,b] : parameters.excess)
//         excess[a]=b;

//     for(auto [v,w,cap, cost] : parameters.EW){
//         edges.push_back({v,w});
//         ep[{v,w}].capacity += cap;
//         ep[{v,w}].cost += cost;
//         ep[{w,v}].cost -= cost;
//     }
//     NetworKit::Graph G = build_graph(parameters.N, edges, false);
//     G.removeMultiEdges();
//     auto algorithm = Koala::SuccessiveApproxMCC(G, ep, excess);
//     algorithm.run();
//     EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
// }

std::tuple<NetworKit::Graph, Koala::edgeid_map<long long>, Koala::node_map<long long>> getInstance(MinCostFlowParams const& params) {
    std::list<std::tuple<int, int, int>> edges;
    Koala::edgeid_map<long long> costs;
    Koala::node_map<long long> b;

    for (auto [key, value] : params.excess) {
        b[key] = value;
    }

    NetworKit::edgeid index = 0;
    for (auto [u, v, capacity, cost] : params.EW) {
        edges.push_back({u, v, capacity});
        costs[index++] = cost;
    }
    NetworKit::Graph G = build_graph(params.N, edges, true, true);

    return { G, costs, b }; 
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

class EdmondsKarpTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(EdmondsKarpTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto [ G, costs, b ] = getInstance(parameters);
    auto algorithm = Koala::EdmondsKarpMCF(G, costs, b);
    algorithm.run();
    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example, EdmondsKarpTest, testing::Values(
    MinCostFlowParams{
        4, {{0, 2, 2, 1}, {2, 0, 1, 0}, {0, 3, 3, 1}, {2, 3, 2, 0}, {1, 2, 2, 1}, {1, 3, 2, 1}}, {{0,3}, {1,2}, {3,-5}}, 5
    },
    MinCostFlowParams{
        4, {{0, 1, 5, 0}, {1, 2, 1, 1}, {1, 2, 2, 2}, {1, 2, 3, 3}, {2, 3, 5, 0}}, {{0, 4}, {3, -4}}, 8  
    },
    MinCostFlowParams{8, {{1, 0, 3, 5}, {2, 0, 2, 1}, {0, 3, 6, 1}, {5, 4, 0, 0}, {6, 4, 4, 3}, {6, 7, 1, 2}, {7, 4, 2, 1}},
        {{7, 2}, {6, 2}, {4, -4}, {0, -1}, {1, 2}, {2, 2}, {3, -3}}, 23
    }, 
    MinCostFlowParams{3, {{0, 1, 5, 1}, {1, 0, 5, 2}, {0, 2, 5, 1}, {2, 0, 5, 2}, {1, 2, 3, 1}, {2, 1, 4, 2}}, 
        {{0, 2}, {1, 3}, {2, -5}}, 5
    }
));

class OrlinTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(OrlinTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto [ G, costs, b ] = getInstance(parameters);
    auto algorithm = Koala::OrlinMCF(G, costs, b);
    algorithm.run();
    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example, OrlinTest, testing::Values(
    MinCostFlowParams{
        4, {{0, 2, 2, 1}, {2, 0, 1, 0}, {0, 3, 3, 1}, {2, 3, 2, 0}, {1, 2, 2, 1}, {1, 3, 2, 1}}, {{0,3}, {1,2}, {3,-5}}, 5
    },
    MinCostFlowParams{
        4, {{0, 1, 5, 0}, {1, 2, 1, 1}, {1, 2, 2, 2}, {1, 2, 3, 3}, {2, 3, 5, 0}}, {{0, 4}, {3, -4}}, 8  
    },
    MinCostFlowParams{8, {{1, 0, 3, 5}, {2, 0, 2, 1}, {0, 3, 6, 1}, {5, 4, 0, 0}, {6, 4, 4, 3}, {6, 7, 1, 2}, {7, 4, 2, 1}},
        {{7, 2}, {6, 2}, {4, -4}, {0, -1}, {1, 2}, {2, 2}, {3, -3}}, 23
    }, 
    MinCostFlowParams{3, {{0, 1, 5, 1}, {1, 0, 5, 2}, {0, 2, 5, 1}, {2, 0, 5, 2}, {1, 2, 3, 1}, {2, 1, 4, 2}}, 
        {{0, 2}, {1, 3}, {2, -5}}, 5
    }
));

class SuccessiveApproxTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(SuccessiveApproxTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto [ G, costs, b ] = getInstance(parameters);
    auto algorithm = Koala::SuccessiveApproxMCC(G, costs, b);
    algorithm.run();
    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example, SuccessiveApproxTest, testing::Values(
    MinCostFlowParams{
        4, {{0, 2, 2, 1}, {2, 0, 1, 0}, {0, 3, 3, 1}, {2, 3, 2, 0}, {1, 2, 2, 1}, {1, 3, 2, 1}}, {{0,3}, {1,2}, {3,-5}}, 5
    },
    MinCostFlowParams{
        4, {{0, 1, 5, 0}, {1, 2, 1, 1}, {1, 2, 2, 2}, {1, 2, 3, 3}, {2, 3, 5, 0}}, {{0, 4}, {3, -4}}, 8  
    },
    MinCostFlowParams{8, {{1, 0, 3, 5}, {2, 0, 2, 1}, {0, 3, 6, 1}, {5, 4, 0, 0}, {6, 4, 4, 3}, {6, 7, 1, 2}, {7, 4, 2, 1}},
        {{7, 2}, {6, 2}, {4, -4}, {0, -1}, {1, 2}, {2, 2}, {3, -3}}, 23
    }, 
    MinCostFlowParams{3, {{0, 1, 5, 1}, {1, 0, 5, 2}, {0, 2, 5, 1}, {2, 0, 5, 2}, {1, 2, 3, 1}, {2, 1, 4, 2}}, 
        {{0, 2}, {1, 3}, {2, -5}}, 5
    }
));