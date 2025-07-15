#include <gtest/gtest.h>

#include <list>
#include <shortest_path/FredericksonPlanarSSSP.hpp>
#include <shortest_path/HenzingerPlanarSSSP.hpp>

#include "helpers.hpp"

struct SSSPParameters {
    std::string name;
    int N;
    std::list<std::tuple<int, int, int>> E;
    NetworKit::node s;
    NetworKit::node t;
    int expected_distance;
    bool directed = false;
};

class PlanarSSSPTest : public testing::TestWithParam<SSSPParameters> {};

TEST_P(PlanarSSSPTest, TestPlanarSSSPSolution) {
    SSSPParameters const& parameters = GetParam();

    std::cout << "\nTest name: " << parameters.name << std::endl;

    NetworKit::Graph G = build_graph(parameters.N, parameters.E, parameters.directed);

    if (parameters.directed) {
        std::cout << "Running HenzingerPlanarSSSP\n" << std::endl;
        Koala::HenzingerPlanarSSSP algorithm(G, parameters.s, parameters.t);
        algorithm.run();
        EXPECT_EQ(algorithm.distance(parameters.t), parameters.expected_distance);
    } else {
        std::cout << "Running FredericksonPlanarSSSP\n" << std::endl;
        Koala::FredericksonPlanarSSSP algorithm(G, parameters.s, parameters.t);
        algorithm.run();
        EXPECT_EQ(algorithm.distance(parameters.t), parameters.expected_distance);
    }
}

SSSPParameters generate_grid_edges(int rows, int cols, bool directed = false) {
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

    return SSSPParameters{"grid", rows * cols, edges, 0,
        static_cast<NetworKit::node>(rows * cols - 1), rows + cols - 2, directed};
}

SSSPParameters generate_line_graph_edges(int length, bool directed = false) {
    std::list<std::tuple<int, int, int>> edges;

    for (int i = 0; i < length - 1; ++i) {
        edges.emplace_back(i, i + 1, 1);
    }

    return {"line", length, edges, 0, length - 1, length - 1, directed};
}

SSSPParameters generate_triangular_grid_edges(int levels, bool directed = false) {
    std::list<std::tuple<int, int, int>> edges;

    int node = 0;
    int row_start;
    int next_row_start;
    int current, left_child, right_child;
    for (int row = 0; row < levels - 1; ++row) {
        row_start = node;
        next_row_start = node + row + 1;

        for (int i = 0; i <= row; ++i) {
            current = row_start + i;
            left_child = next_row_start + i;
            right_child = next_row_start + i + 1;

            edges.emplace_back(current, left_child, 1);
            edges.emplace_back(current, right_child, 1);
            edges.emplace_back(left_child, right_child, 1);
        }

        node = next_row_start;
    }

    return {"triangular-grid", right_child + 1, edges, 0, right_child, levels - 1, directed};
}

INSTANTIATE_TEST_SUITE_P(Default, PlanarSSSPTest,
    testing::Values(generate_grid_edges(40, 40), generate_line_graph_edges(100),
        generate_triangular_grid_edges(10), generate_grid_edges(40, 40, true),
        generate_line_graph_edges(100, true), generate_triangular_grid_edges(10, true)));
