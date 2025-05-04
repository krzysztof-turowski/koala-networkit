#include <gtest/gtest.h>
#include <list>

#include "planar_sssp/PlanarSSSP.hpp"
#include "planar_sssp/FredericksonPlanarSSSP.hpp"
#include "helpers.hpp"

struct SSSPParameters {
    std::string name;
    int N;
    std::list<std::tuple<int, int, int>> E;
    NetworKit::node s;
    NetworKit::node t;
    int expectedDistance;
};

class PlanarSSSPTest : public testing::TestWithParam<SSSPParameters> {
};

TEST_P(PlanarSSSPTest, TestPlanarSSSPSolution) {
    SSSPParameters const& parameters = GetParam();

    std::cout << "Test name: " << parameters.name << std::endl;

    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);

    Koala::FredericksonPlanarSSSP algorithm(G, parameters.s, parameters.t);

    algorithm.run();

    EXPECT_EQ(algorithm.getSSSDistanceToTarget(), parameters.expectedDistance);
}

std::list<std::tuple<int, int, int>> generateGridEdges(int rows, int cols) {
    std::list<std::tuple<int, int, int>> edges;

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int node = r * cols + c;
            if (c < cols - 1) {
                edges.emplace_back(node, node + 1, 1);
            }
            if (r < rows - 1) {
                edges.emplace_back(node, node + cols, 1);
            }
        }
    }

    return edges;
}

std::list<std::tuple<int, int, int>> generateLineGraphEdges(int length) {
    std::list<std::tuple<int, int, int>> edges;

    for (int i = 0; i < length - 1; ++i) {
        edges.emplace_back(i, i + 1, 1);
    }

    return edges;
}

std::list<std::tuple<int, int, int>> generateTriangularGridEdges(int levels) {
    std::list<std::tuple<int, int, int>> edges;

    int node = 0; // numeracja węzłów
    for (int row = 0; row < levels - 1; ++row) {
        int rowStart = node;
        int nextRowStart = node + row + 1;

        for (int i = 0; i <= row; ++i) {
            int current = rowStart + i;
            int leftChild = nextRowStart + i;
            int rightChild = nextRowStart + i + 1;

            edges.emplace_back(current, leftChild, 1);
            edges.emplace_back(current, rightChild, 1);
            edges.emplace_back(leftChild, rightChild, 1);
        }

        node = nextRowStart;
    }

    return edges;
}



INSTANTIATE_TEST_SUITE_P(
    Default, PlanarSSSPTest,testing::Values(
    // SSSPParameters{ "star", 5, {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 6}}, 0, 2, 5 },
    // SSSPParameters{ "double-star", 8, {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 6}, {4,7,1}, {4,6,1}, { 4,5,1 }}, 0, 2, 5 },
    SSSPParameters{ "10x10-grid", 625, generateGridEdges(25, 25), 0, 99, 18 },
    SSSPParameters{ "line", 1000, generateLineGraphEdges(1000), 0, 99, 99 },
    SSSPParameters{ "triangular-grid", 55, generateTriangularGridEdges(10), 0, 45, 18 }

        // SSSPParameters{ "", 4, {{0, 1, 1}, {0, 3, 2}, {1, 2, 3}, {2, 3, 4}}, 0, 2, 5 },

        // SSSPParameters{ "", 4, {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {1, 2, 1}, {1, 3, 1}, {2, 3, 1}}, 0, 2, 5 },

        // SSSPParameters{ "", 5, {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {0, 4, 1}}, 0, 2, 5 },

        // SSSPParameters{ "", 4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 3}}, 0, 2, 5 },

        // SSSPParameters{ "", 5, {{0, 1, 10}, {0, 4, 3}, {1, 2, 2}, {2, 3, 4}, {3, 4, 5}}, 0, 2, 5 }
    ));
    