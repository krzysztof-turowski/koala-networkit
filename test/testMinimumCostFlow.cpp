// TODO: undirected graph, antisymmetric costs
#include <gtest/gtest.h>

#include <list>

#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>
#include <flow/minimum_cost_flow/EdmondsKarpMCF.hpp>
#include <flow/minimum_cost_flow/OrlinMCF.hpp>

#include "helpers.hpp"

// #define DEBUG_DUMP
#ifdef DEBUG_DUMP
#define DBG(x) std::cerr << x
#else 
#define DBG(x)
#endif

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

Koala::MCFlowNetwork getInstance(MinCostFlowParams const& params) {
    std::list<std::tuple<int, int, int>> edges;

    std::unordered_map<NetworKit::Edge, std::int64_t> costs;
    std::unordered_map<NetworKit::node, std::int64_t> excess(params.excess.begin(), params.excess.end());

    for (auto [u, v, capacity, cost] : params.EW) {
        costs[{u,v}] = cost;
        edges.push_back({u, v, capacity});
    }
    NetworKit::Graph G = build_graph(params.N, edges, true, false);

    return Koala::MCFlowNetwork(G, costs, excess); 
}

const std::vector<MinCostFlowParams> basic_tests = {
    MinCostFlowParams{
        4, {{0, 2, 2, 1}, {2, 0, 1, 0}, {0, 3, 3, 1}, {2, 3, 2, 0}, {1, 2, 2, 1}, {1, 3, 2, 1}}, {{0,3}, {1,2}, {3,-5}}, 5
    },
    MinCostFlowParams{8, {{1, 0, 3, 5}, {2, 0, 2, 1}, {0, 3, 6, 1}, {5, 4, 0, 0}, {6, 4, 4, 3}, {6, 7, 1, 2}, {7, 4, 2, 1}},
        {{7, 2}, {6, 2}, {4, -4}, {0, -1}, {1, 2}, {2, 2}, {3, -3}}, 23
    }, 
    MinCostFlowParams{
        3, 
        {{0, 1, 5, 1}, {1, 2, 5, 2}}, 
        {{0, 3}, {2, -3}}, 
        9
    },
    MinCostFlowParams{
        4, 
        {{0, 1, 2, 1}, {1, 3, 2, 1}, {0, 2, 5, 5}, {2, 3, 5, 5}}, 
        {{0, 4}, {3, -4}}, 
        24
    },
    MinCostFlowParams{
        2, 
        {{0, 1, 10, 3}}, 
        {{0, 2}, {1, -2}}, 
        6
    },
    MinCostFlowParams{
        5, 
        {{0, 2, 3, 2}, {1, 2, 4, 1}, {2, 3, 5, 3}, {2, 4, 5, 4}}, 
        {{0, 2}, {1, 3}, {3, -1}, {4, -4}}, 
        26
    },
    MinCostFlowParams{
        4, 
        {{0, 1, 5, 0}, {1, 2, 5, 0}, {2, 3, 5, 0}, {0, 3, 1, 10}}, 
        {{0, 4}, {3, -4}}, 
        0
    },
    MinCostFlowParams{
        4,
        {{0, 1, 10, 1}, {1, 2, 5, 5}, {2, 3, 10, 1}},
        {{0, 5}, {3, -5}},
        35
    }
};

const std::vector<MinCostFlowParams> complex_tests = {
    MinCostFlowParams{
        6, 
        {
            {0, 3, 1, 10}, {0, 4, 1, 2}, {0, 5, 1, 8},
            {1, 3, 1, 9},  {1, 4, 1, 8}, {1, 5, 1, 1},
            {2, 3, 1, 2},  {2, 4, 1, 9}, {2, 5, 1, 8}
        }, 
        {{0, 1}, {1, 1}, {2, 1}, {3, -1}, {4, -1}, {5, -1}}, 
        5
    },
    MinCostFlowParams{
        4,
        {
            {0, 1, 1, 0}, {0, 2, 1, 2}, 
            {1, 2, 1, 1}, {1, 3, 1, 4}, 
            {2, 3, 1, 0}
        },
        {{0, 2}, {3, -2}},
        6
    },
    MinCostFlowParams{
        6,
        {
            {0, 1, 10, 2}, {0, 2, 5, 5},
            {1, 2, 5, 1},  {1, 3, 5, 4},
            {2, 3, 3, 1},  {2, 4, 8, 2},
            {3, 4, 5, 1},  {3, 5, 10, 3},
            {4, 5, 10, 2}
        },
        {{0, 12}, {5, -12}},
        98
    },
    MinCostFlowParams{
        5,
        {
            {0, 1, 10, 100}, {0, 2, 10, 1},
            {2, 3, 10, 1},   {3, 1, 10, 1},
            {1, 4, 10, 1}
        },
        {{0, 5}, {4, -5}},
        20
    },
    MinCostFlowParams{
        7,
        {
            {0, 3, 15, 2}, {1, 3, 10, 3}, {2, 3, 10, 1},
            {0, 4, 5, 10},
            {3, 4, 10, 2}, {3, 5, 10, 4}, {3, 6, 10, 3},
            {1, 5, 5, 6},  {2, 6, 5, 2}
        },
        {{0, 10}, {1, 5}, {2, 5}, {4, -8}, {5, -7}, {6, -5}},
        84
    }
};

/*
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
*/

class EdmondsKarpTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(EdmondsKarpTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto network = getInstance(parameters);
    auto algorithm = Koala::EdmondsKarpMCF(network);
    algorithm.run();

    network.getGraph().forEdges([&](NetworKit::node u, NetworKit::node v) {
        DBG(u << ' ' << v << " -> " << algorithm.getFlow({u, v}) << '\n');
    });

    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example_edmonds, EdmondsKarpTest, testing::ValuesIn(basic_tests));
INSTANTIATE_TEST_SUITE_P(test_complex_edmonds, EdmondsKarpTest, testing::ValuesIn(complex_tests));

class OrlinTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(OrlinTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto network = getInstance(parameters);
    auto algorithm = Koala::OrlinMCF(network);
    algorithm.run();

    network.getGraph().forEdges([&](NetworKit::node u, NetworKit::node v) {
        DBG(u << ' ' << v << " -> " << algorithm.getFlow({u, v}) << '\n');
    });

    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example_orlin, OrlinTest, testing::ValuesIn(basic_tests));
INSTANTIATE_TEST_CASE_P(test_complex_orlin, OrlinTest, testing::ValuesIn(complex_tests));



class SuccessiveApproxTest
    : public testing::TestWithParam<MinCostFlowParams> { };

TEST_P(SuccessiveApproxTest, test) {
    MinCostFlowParams const& parameters = GetParam();
    auto network = getInstance(parameters);
    auto algorithm = Koala::SuccessiveApproxMCC(network);
    algorithm.run();

    network.getGraph().forEdges([&](NetworKit::node u, NetworKit::node v) {
        DBG(u << ' ' << v << " -> " << algorithm.getFlow({u, v}) << '\n');
    });

    EXPECT_EQ(algorithm.getMinCost(), parameters.minCost);
}

INSTANTIATE_TEST_SUITE_P(test_example_successive, SuccessiveApproxTest, testing::ValuesIn(basic_tests));
INSTANTIATE_TEST_SUITE_P(test_complex_successive, SuccessiveApproxTest, testing::ValuesIn(complex_tests));