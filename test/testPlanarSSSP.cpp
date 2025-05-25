#include <gtest/gtest.h>
#include <list>

#include "planar_sssp/PlanarSSSP.hpp"
#include "planar_sssp/FredericksonPlanarSSSP.hpp"
#include "planar_sssp/HenzingerPlanarSSSP.hpp"
#include "helpers.hpp"

struct SSSPParameters {
    std::string name;
    int N;
    std::list<std::tuple<int, int, int>> E;
    NetworKit::node s;
    NetworKit::node t;
    int expectedDistance;
    bool directed = false;
};

class PlanarSSSPTest : public testing::TestWithParam<SSSPParameters> {
};

TEST_P(PlanarSSSPTest, TestPlanarSSSPSolution) {
    SSSPParameters const& parameters = GetParam();

    std::cout << "\nTest name: " << parameters.name << std::endl;

    NetworKit::Graph G = build_graph(parameters.N, parameters.E, parameters.directed);

    if(parameters.directed){
        std::cout << "Running HenzingerPlanarSSSP\n" << std::endl;
        Koala::HenzingerPlanarSSSP algorithm(G, parameters.s, parameters.t);
        algorithm.run();
        EXPECT_EQ(algorithm.getSSSDistanceToTarget(), parameters.expectedDistance);
    }else{
        std::cout<< "Running FredericksonPlanarSSSP\n" << std::endl;
        Koala::FredericksonPlanarSSSP algorithm(G, parameters.s, parameters.t);
        algorithm.run();
        EXPECT_EQ(algorithm.getSSSDistanceToTarget(), parameters.expectedDistance);
    }
}

SSSPParameters generateGridEdges(int rows, int cols, bool directed = false) {
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

    return SSSPParameters{ "grid", rows*cols, edges, 0, static_cast<NetworKit::node>(rows*cols-1), rows + cols - 2, directed };
}

SSSPParameters generateLineGraphEdges(int length, bool directed = false) {
    std::list<std::tuple<int, int, int>> edges;

    for (int i = 0; i < length - 1; ++i) {
        edges.emplace_back(i, i + 1, 1);
    }

    return {"line", length, edges, 0, length-1, length-1, directed};
}

SSSPParameters generateTriangularGridEdges(int levels, bool directed = false) {
    std::list<std::tuple<int, int, int>> edges;

    int node = 0; // numeracja węzłów
    int rowStart;
    int nextRowStart;
    int current;
    int leftChild, rightChild; 
    for (int row = 0; row < levels - 1; ++row) {
        rowStart = node;
        nextRowStart = node + row + 1;

        for (int i = 0; i <= row; ++i) {
            current = rowStart + i;
            leftChild = nextRowStart + i;
            rightChild = nextRowStart + i + 1;

            edges.emplace_back(current, leftChild, 1);
            edges.emplace_back(current, rightChild, 1);
            edges.emplace_back(leftChild, rightChild, 1);
        }

        node = nextRowStart;
    }

    return {"triangular-grid", rightChild + 1, edges, 0, rightChild, levels-1, directed};

}

INSTANTIATE_TEST_SUITE_P(
    Default, PlanarSSSPTest,testing::Values(
    generateGridEdges(40, 40),
    generateLineGraphEdges(100),
    generateTriangularGridEdges(10),
    generateGridEdges(40, 40, true),
    generateLineGraphEdges(100, true),
    generateTriangularGridEdges(10, true)
    ));
    