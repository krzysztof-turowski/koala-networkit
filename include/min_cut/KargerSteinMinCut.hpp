/*
 * KargerSteinMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using the Karger-Stein algorithm.
 * Created on: 05.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Min-Cut problem on a given graph
 * using the Karger-Stein algorithm.
 */
class KargerSteinMinCut {
public:
    /**
     * Constructor that initializes the graph from a given adjacency matrix.
     *
     * @param vertices The number of vertices in the graph.
     * @param repeats The number of repetitions to perform the algorithm.
     * @param graphMatrix The graph represented as an adjacency matrix.
     */
    KargerSteinMinCut(int vertices, int repeats, const std::vector<std::vector<int>>& graphMatrix);

    /**
     * Adds an edge to the graph array.
     *
     * @param u The starting vertex of the edge.
     * @param v The ending vertex of the edge.
     * @param weight The weight of the edge.
     */
    void addEdge(int u, int v, int weight);

    /**
     * Executes the Min-Cut problem solver using the Karger-Stein algorithm.
     */
    void solve();

    /**
     * Retrieves the minimum cut value found by the algorithm.
     *
     * @return The minimum cut value.
     */
    int getMinCut();

private:
    int numVertices, numRepeats;
    std::vector<NetworKit::Edge> edges;
    int bestMinCut;

    int find(std::vector<int>& parent, int i);
    void unionSub(std::vector<int>& parent, std::vector<int>& rank, int x, int y);
    int recursiveMinCut(int numVertices);
    int standardMinCut();
};

} // namespace Koala