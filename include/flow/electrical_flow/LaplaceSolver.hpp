#pragma once

#include <eigen3/Eigen/Core>

#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * Solve a weighted graph Laplacian system.
 *
 * @param graph Input graph.
 * @param weights Conductance matrix indexed by edge endpoints.
 * @param b Demand vector.
 * @return A potential vector satisfying the Laplacian system in the pseudoinverse sense.
 */
Eigen::VectorXd solve_laplace(
    const NetworKit::Graph &graph, const std::vector<std::vector<double>> &weights,
    const Eigen::VectorXd &b);

}  // namespace Koala
