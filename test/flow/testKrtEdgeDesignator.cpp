#include <gtest/gtest.h>

#include <stdexcept>

#include "flow/maximum_flow/KrtEdgeDesignator.hpp"

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
