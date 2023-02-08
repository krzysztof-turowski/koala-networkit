#include <gtest/gtest.h>

#include <list>

#include <coloring/EnumerationVertexColoring.hpp>

struct VertexColoringParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int colors;
};

class BrownsOrdinaryEnumerationVertexColoringTest
    : public testing::TestWithParam<VertexColoringParameters> { };

NetworKit::Graph build_graph(const int& N, const std::list<std::pair<int, int>>& E) {
    NetworKit::Graph G(N, false, false);
    for (const auto& [u, v] : E) {
        G.addEdge(u, v);
    }
    return G;
}

TEST_P(BrownsOrdinaryEnumerationVertexColoringTest, test) {
    VertexColoringParameters const& parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E);
    auto algorithm = Koala::BrownsOrdinaryEnumerationVertexColoring(G);
    algorithm.run();

    auto colors = algorithm.getColoring();

    for (const auto& [u, v] : parameters.E) {
        EXPECT_NE(colors[u], colors[v]);
    }

    int max_color = 0;
    for (const auto& [v, c] : colors) {
        max_color = std::max(max_color, c);
    }
    EXPECT_EQ(max_color, parameters.colors);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, BrownsOrdinaryEnumerationVertexColoringTest, testing::Values(
        VertexColoringParameters{ 4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, 2 },
        VertexColoringParameters{ 6, {{0, 1}, {0, 2}, {1, 2}, {0, 3}, {3, 4}, {1, 5}, {4, 5}}, 3 },
        VertexColoringParameters{
            10, {{0, 2}, {0, 3}, {0, 4}, {0, 5}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 6}, {2, 8},
                 {3, 4}, {3, 6}, {4, 6}, {5, 6}, {5, 7}, {5, 9}, {7, 8}, {7, 9}, {8, 9}}, 4},
        VertexColoringParameters{
            8, {{0, 1}, {0, 2}, {0, 4}, {0, 6}, {1, 2}, {1, 3}, {1, 7}, {2, 3}, {2, 4}, {3, 5},
                {3, 7}, {4, 5}, {4, 6}, {5, 6}, {5, 7}, {6, 7}}, 4},
        VertexColoringParameters{
            8, {{0, 1}, {0, 2}, {0, 4}, {0, 6}, {1, 3}, {1, 5}, {1, 7}, {2, 3}, {2, 4}, {2, 5},
                {2, 6}, {3, 4}, {3, 5}, {3, 7}, {4, 6}, {4, 7}, {5, 6}, {5, 7}, {6, 7}}, 5},
        VertexColoringParameters{
            10, {{0, 1}, {0, 2}, {0, 3}, {0, 5}, {0, 6}, {0, 7}, {1, 2}, {1, 3}, {1, 4}, {1, 6},
                 {1, 7}, {2, 3}, {2, 4}, {2, 5}, {2, 7}, {3, 4}, {3, 5}, {3, 6}, {4, 5}, {4, 7},
                 {4, 8}, {4, 9}, {5, 6}, {5, 8}, {5, 9}, {6, 7}, {6, 8}, {6, 9}, {7, 8}, {7, 9},
                 {8, 9}}, 5}
));
