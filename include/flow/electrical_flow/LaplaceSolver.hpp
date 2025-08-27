#pragma once

#include <eigen3/Eigen/Core>
#include <networkit/graph/Graph.hpp>

using namespace std;
using namespace Eigen;
using namespace NetworKit;

namespace Koala {
VectorXd solveLaplace(const Graph &graph, const vector<vector<double>> &weights,
                      const VectorXd &b);
}
