#pragma once

#include <algorithm>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup shortest_path
 * The base class for the All-Pairs Shortest Paths algorithms.
 *
 * Holds the pairwise distance matrix and exposes the getters shared by every
 * APSP algorithm, together with the diameter (the largest finite distance) --
 * a computation common to all of them. Concrete algorithms fill in
 * distances inside run() and then call computeDiameter().
 *
 * @tparam WeightType type storing a single distance (e.g. @c int for
 *         unweighted graphs, @c NetworKit::edgeweight for weighted graphs).
 *         Unreachable pairs are stored as std::numeric_limits<WeightType>::max().
 */
template <typename WeightType>
class AllPairsShortestPaths : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the all-pairs shortest paths procedure.
     *
     * @param graph The input graph.
     */
    explicit AllPairsShortestPaths(const NetworKit::Graph &graph) : graph(graph) { }

    /**
     * Given an input graph to consume, set up the all-pairs shortest paths
     * procedure without copying it.
     *
     * @param graph The input graph, moved into the algorithm.
     */
    explicit AllPairsShortestPaths(NetworKit::Graph &&graph) : graph(std::move(graph)) { }

    void run() override = 0;

    /**
     * Return the distance matrix found by the algorithm.
     *
     * @return distances[u][v] is the estimated distance between u and v.
     */
    const std::vector<std::vector<WeightType>>& getDistances() const {
        assureFinished();
        return distances;
    }

    /**
     * Return the distance between a single pair of nodes.
     */
    WeightType getDistance(NetworKit::node u, NetworKit::node v) const {
        assureFinished();
        return distances[u][v];
    }

    /**
     * Return the diameter: the largest finite pairwise distance.
     */
    WeightType getDiameter() const {
        assureFinished();
        return diameter;
    }

 protected:
    /**
     * Derive the diameter from the distance matrix, ignoring unreachable pairs
     * (distance == std::numeric_limits<WeightType>::max()). Call once the
     * matrix has been filled in run().
     */
    void computeDiameter() {
        constexpr WeightType INF = std::numeric_limits<WeightType>::max();
        diameter = WeightType{};
        for (const auto &row : distances) {
            for (const WeightType d : row) {
                if (d < INF) {
                    diameter = std::max(diameter, d);
                }
            }
        }
    }

    std::optional<NetworKit::Graph> graph;
    std::vector<std::vector<WeightType>> distances;
    WeightType diameter = WeightType{};
};

}  // namespace Koala
