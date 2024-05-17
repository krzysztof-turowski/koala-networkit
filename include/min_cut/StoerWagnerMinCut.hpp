/*
 * StoerWagnerMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using the Stoer-Wagner algorithm.
 * Created on: 07.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Min-Cut problem on a given graph
 * using the Stoer-Wagner algorithm.
 */
class StoerWagnerMinCut {
public:
    /**
     * Constructor that initializes the graph from a given adjacency matrix.
     *
     * @param vertices The number of vertices in the graph.
     * @param graphMatrix The graph represented as an adjacency matrix.
     */
    StoerWagnerMinCut(int vertices, const std::vector<std::vector<int>>& graphMatrix);

    /**
     * Executes the Min-Cut problem solver using the Stoer-Wagner algorithm.
     */
    void solve();

    /**
     * Retrieves the minimum cut value found by the algorithm.
     *
     * @return The minimum cut value.
     */
    int getMinCut() const;

private:
    int numVertices;
    std::vector<std::vector<int>> graph;
    std::vector<int> vertices;
    int bestMinCut;

    int minCutPhase(std::vector<int>& vertices, std::vector<std::vector<int>>& weights);
    int findMinCut(std::vector<std::vector<int>>& weights);
};

} // namespace Koala