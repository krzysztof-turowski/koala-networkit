/*
 * HaoOrlinMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using the Hao-Orlin algorithm.
 * Created on: 30.04.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Min-Cut problem on a given graph
 * using the Hao-Orlin algorithm.
 */
class HaoOrlinMinCut {
public:
    /**
     * Constructor that initializes the solver with a graph.
     *
     * @param vertices The number of vertices in the graph.
     * @param src The source vertex.
     * @param snk The sink vertex.
     * @param graphMatrix The graph represented as an adjacency matrix.
     */
    HaoOrlinMinCut(int vertices, int src, int snk, const std::vector<std::vector<int>>& graphMatrix);

    /**
     * Executes the Min-Cut problem solver using the Hao-Orlin algorithm.
     */
    void findMinimumCut();

    /**
     * Retrieves the minimum cut value found by the algorithm.
     *
     * @return The minimum cut value.
     */
    long long getMinCut() const;

private:
    int numVertices;
    int source;
    int sink;
    std::vector<std::vector<int>> graph;
    int bestValue;
    std::vector<int> cut;
};

} // namespace Koala
