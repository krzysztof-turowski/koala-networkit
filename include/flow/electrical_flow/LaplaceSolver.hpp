#pragma once

#include <Eigen/Core>
#include <networkit/graph/Graph.hpp>

namespace Koala {
Eigen::VectorXd solveLaplace(const NetworKit::Graph &graph,
                             const std::vector<std::vector<double>> &weights,
                             const Eigen::VectorXd &b);
}
