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
 * @ingroup min-cut
 * The class for solving the Min-Cut problem on a given graph
 * using the Push-Relabel algorithm.
 */
class PushRelabelMinCut {
 public:
    using MinCut::MinCut;

    /**
     * Executes the Min-Cut problem solver.
     */
    void run();

 private:
    struct Edge {
        int source, sink;
        int capacity, flow;
        int revIndex;
    };

    int numVertices;
    int source;
    int sink;
    std::vector<std::vector<Edge>> graph;
    std::vector<int> excess;              // Excess flow at each vertex
    std::vector<int> height;              // Height of each vertex for the push-relabel method
    std::vector<bool> active;             // Active vertices for optimization
    std::vector<int> current;             // Current edge being considered for each vertex

    // Helper methods
    void initializePreflow();
    void push(Edge& e);
    void relabel(int v);
    void discharge(int v);
    void addEdge(int from, int to, int capacity);
    int overflowVertex();
};

}  // namespace Koala
