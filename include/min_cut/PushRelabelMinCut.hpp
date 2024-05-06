/*
 * PushRelabelMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using the Push-Relabel algorithm.
 * Created on: 30.04.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Min-Cut problem on a given graph
 * using the Push-Relabel algorithm.
 */
class PushRelabelMinCut {
public:
    /**
     * Constructor that initializes the solver with a graph.
     *
     * @param vertices The number of vertices in the graph.
     * @param src The source vertex.
     * @param snk The sink vertex.
     * @param graphMatrix The graph represented as an adjacency matrix.
     */
    PushRelabelMinCut(int vertices, int src, int snk, const std::vector<std::vector<int>>& graphMatrix);

    /**
     * Executes the Min-Cut problem solver using the Push-Relabel algorithm.
     */
    void solve();

    /**
     * Retrieves the maximum flow value found by the algorithm.
     *
     * @return The maximum flow value.
     */
    long long getMaxFlow() const;

    /**
     * Retrieves the best set partition found by the algorithm.
     * 1 represents vertices in one set, while 0 represents vertices in the opposite set.
     *
     * @return Vector of integers indicating the set membership of each vertex.
     */
    std::vector<int> getMinCut() const;

private:
    struct Edge {
        int source, sink;
        long long capacity, flow;
        int revIndex;
    };

    int numVertices;
    int source;
    int sink;
    std::vector<std::vector<Edge>> graph;
    std::vector<long long> excess;        // Excess flow at each vertex
    std::vector<int> height;              // Height of each vertex for the push-relabel method
    std::vector<bool> active;             // Active vertices for optimization
    std::vector<int> current;             // Current edge being considered for each vertex

    // Helper methods
    void initializePreflow();
    void push(Edge& e);
    void relabel(int v);
    void discharge(int v);
    void addEdge(int from, int to, long long capacity);
    int overflowVertex();
};

} // namespace Koala
