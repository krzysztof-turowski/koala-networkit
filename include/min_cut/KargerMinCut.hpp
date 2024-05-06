/*
 * KargerMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using Karger's algorithm.
 * Created on: 03.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Min-Cut problem on a given graph
 * using Karger's algorithm.
 */
class KargerMinCut {
public:
    /**
     * Constructor that initializes the graph from a given adjacency matrix.
     *
     * @param vertices The number of vertices in the graph.
     * @param repeat The number of repetitions to perform the algorithm.
     * @param graphMatrix The graph represented as an adjacency matrix.
     */
    KargerMinCut(int vertices, int repeat, const std::vector<std::vector<int>>& graphMatrix);

    /**
     * Adds an edge to the graph array.
     *
     * @param u The starting vertex of the edge.
     * @param v The ending vertex of the edge.
     * @param weight The weight of the edge.
     */
    void addEdge(int u, int v, int weight);

    /**
     * Executes the Min-Cut problem solver using Karger's algorithm.
     */
    void solve();

    /**
     * Retrieves the minimum cut value found by the algorithm.
     *
     * @return The minimum cut value.
     */
    int getMinCut();

private:
    int numVertices, numRepeat;
    std::vector<NetworKit::Edge> edges;
    int bestMinCut;

    int findMinCut();
    int find(std::vector<int>& parent, int i);
    void unionSub(std::vector<int>& parent, std::vector<int>& rank, int x, int y);
};

} // namespace Koala
